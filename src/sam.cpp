#include "burnham/sam.hpp"

#include "burnham/cigar.hpp"
#include "burnham/fs.hpp"
#include "burnham/text.hpp"

#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <map>
#include <optional>
#include <sstream>
#include <string_view>
#include <vector>

namespace burnham {
namespace {

struct SamRecordLine {
    std::string line;
    std::vector<std::string> fields;
    std::int64_t pos = 0;
};

struct OrderedContig {
    std::string name;
    std::int64_t length = 0;
};

struct DuplicateEndpoint {
    std::string contig;
    std::int64_t position = 0;
    bool reverse = false;
};

struct DuplicateTemplate {
    std::vector<std::size_t> record_indices;
    std::string key;
    std::int64_t base_quality_score = 0;
    std::int64_t mapq_score = 0;
    std::size_t input_order = 0;
};

bool is_unmapped_flag(const std::string& flag_text) {
    std::int64_t flag = 0;
    return parse_i64(flag_text, flag) && ((flag & 0x4) != 0);
}

std::int64_t parse_flag_or_zero(const std::string& flag_text) {
    std::int64_t flag = 0;
    parse_i64(flag_text, flag);
    return flag;
}

std::string metrics_json(int total, int lifted, int rejected) {
    std::ostringstream out;
    out << "{\n";
    out << "  \"schema_version\": 1,\n";
    out << "  \"records\": " << total << ",\n";
    out << "  \"lifted\": " << lifted << ",\n";
    out << "  \"missing_mapping\": " << rejected << "\n";
    out << "}\n";
    return out.str();
}

std::string read_group_header(const ReadGroupOptions& read_group) {
    std::ostringstream out;
    out << "@RG\tID:" << read_group.id << "\tSM:" << read_group.sample;
    if (!read_group.library.empty()) {
        out << "\tLB:" << read_group.library;
    }
    if (!read_group.platform.empty()) {
        out << "\tPL:" << read_group.platform;
    }
    if (!read_group.platform_unit.empty()) {
        out << "\tPU:" << read_group.platform_unit;
    }
    return out.str();
}

void remove_read_group_tag(std::vector<std::string>& fields) {
    fields.erase(std::remove_if(fields.begin() + std::min<std::size_t>(fields.size(), 11), fields.end(), [](const std::string& field) {
                     return starts_with(field, "RG:Z:");
                 }),
                 fields.end());
}

void remove_tags_with_prefixes(std::vector<std::string>& fields, std::initializer_list<std::string_view> prefixes) {
    fields.erase(std::remove_if(fields.begin() + std::min<std::size_t>(fields.size(), 11), fields.end(), [&](const std::string& field) {
                     for (const auto prefix : prefixes) {
                         if (starts_with(field, prefix)) {
                             return true;
                         }
                     }
                     return false;
                 }),
                 fields.end());
}

void mark_record_unmapped(std::vector<std::string>& fields, const std::string& reason) {
    auto flag = parse_flag_or_zero(fields[1]);
    flag |= 0x4;
    fields[1] = std::to_string(flag);
    fields[2] = "*";
    fields[3] = "0";
    fields[4] = "0";
    fields[5] = "*";
    fields[6] = "*";
    fields[7] = "0";
    fields[8] = "0";
    remove_tags_with_prefixes(fields, {"BC:Z:"});
    fields.push_back("BC:Z:" + reason);
}

Result<std::vector<OrderedContig>> read_ordered_fai(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        return make_error("failed to open FAI: " + path.string());
    }
    std::vector<OrderedContig> contigs;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 2) {
            return make_error("invalid FAI line " + std::to_string(line_number));
        }
        std::int64_t length = 0;
        if (!parse_i64(fields[1], length) || length < 0) {
            return make_error("invalid FAI length on line " + std::to_string(line_number));
        }
        contigs.push_back(OrderedContig{fields[0], length});
    }
    return contigs;
}

std::optional<std::int64_t> record_reference_end(const std::vector<std::string>& fields) {
    if (fields.size() < 6 || fields[3] == "0" || fields[5] == "*") {
        return std::nullopt;
    }
    std::int64_t pos = 0;
    if (!parse_i64(fields[3], pos) || pos <= 0) {
        return std::nullopt;
    }
    auto cigar = parse_cigar(fields[5]);
    if (!cigar) {
        return std::nullopt;
    }
    return pos + cigar_stats(cigar.value()).reference_span - 1;
}

bool record_is_mapped(const std::vector<std::string>& fields) {
    return fields.size() >= 11 && (parse_flag_or_zero(fields[1]) & 0x4) == 0 && fields[2] != "*";
}

void set_or_clear_flag(std::int64_t& flag, std::int64_t mask, bool value) {
    if (value) {
        flag |= mask;
    } else {
        flag &= ~mask;
    }
}

void add_mate_tags(std::vector<std::string>& record, const std::vector<std::string>& mate, bool mate_mapped) {
    remove_tags_with_prefixes(record, {"MC:Z:", "MQ:i:"});
    if (!mate_mapped) {
        return;
    }
    record.push_back("MC:Z:" + mate[5]);
    record.push_back("MQ:i:" + mate[4]);
}

std::int64_t clipped_bases_at_start(const std::vector<CigarOp>& cigar) {
    std::int64_t clipped = 0;
    for (const auto& op : cigar) {
        if (op.op == 'S' || op.op == 'H') {
            clipped += op.length;
            continue;
        }
        break;
    }
    return clipped;
}

std::int64_t clipped_bases_at_end(const std::vector<CigarOp>& cigar) {
    std::int64_t clipped = 0;
    for (auto iter = cigar.rbegin(); iter != cigar.rend(); ++iter) {
        if (iter->op == 'S' || iter->op == 'H') {
            clipped += iter->length;
            continue;
        }
        break;
    }
    return clipped;
}

std::optional<DuplicateEndpoint> duplicate_endpoint(const std::vector<std::string>& fields) {
    if (fields.size() < 11 || !record_is_mapped(fields)) {
        return std::nullopt;
    }
    std::int64_t pos = 0;
    if (!parse_i64(fields[3], pos) || pos <= 0) {
        return std::nullopt;
    }
    auto cigar = parse_cigar(fields[5]);
    if (!cigar || cigar.value().empty()) {
        return std::nullopt;
    }
    const auto flag = parse_flag_or_zero(fields[1]);
    const bool reverse = (flag & 0x10) != 0;
    DuplicateEndpoint endpoint;
    endpoint.contig = fields[2];
    endpoint.reverse = reverse;
    if (reverse) {
        endpoint.position = pos + cigar_stats(cigar.value()).reference_span + clipped_bases_at_end(cigar.value()) - 1;
    } else {
        endpoint.position = pos - clipped_bases_at_start(cigar.value());
    }
    return endpoint;
}

std::string endpoint_key(const DuplicateEndpoint& endpoint) {
    std::ostringstream out;
    out << endpoint.contig << ':' << endpoint.position << ':' << (endpoint.reverse ? '-' : '+');
    return out.str();
}

std::int64_t base_quality_score(const std::vector<std::string>& fields) {
    if (fields.size() < 11 || fields[10] == "*") {
        return 0;
    }
    std::int64_t score = 0;
    for (char quality_char : fields[10]) {
        const auto quality = static_cast<int>(quality_char) - 33;
        if (quality >= 15) {
            score += quality;
        }
    }
    return score;
}

std::int64_t mapq_score(const std::vector<std::string>& fields) {
    if (fields.size() < 5) {
        return 0;
    }
    std::int64_t mapq = 0;
    parse_i64(fields[4], mapq);
    return mapq;
}

bool duplicate_record_eligible(const std::vector<std::string>& fields) {
    if (fields.size() < 11 || !record_is_mapped(fields)) {
        return false;
    }
    const auto flag = parse_flag_or_zero(fields[1]);
    if ((flag & 0x100) != 0 || (flag & 0x800) != 0) {
        return false;
    }
    auto cigar = parse_cigar(fields[5]);
    return cigar.ok() && !cigar.value().empty();
}

std::map<std::string, std::string> read_group_libraries(const std::vector<std::string>& headers) {
    std::map<std::string, std::string> libraries;
    for (const auto& header : headers) {
        if (!starts_with(header, "@RG")) {
            continue;
        }
        auto fields = split(header, '\t');
        std::string id;
        std::string library;
        for (const auto& field : fields) {
            if (starts_with(field, "ID:")) {
                id = field.substr(3);
            } else if (starts_with(field, "LB:")) {
                library = field.substr(3);
            }
        }
        if (!id.empty()) {
            libraries[id] = library.empty() ? id : library;
        }
    }
    return libraries;
}

std::string record_read_group(const std::vector<std::string>& fields) {
    for (std::size_t i = 11; i < fields.size(); ++i) {
        if (starts_with(fields[i], "RG:Z:")) {
            return fields[i].substr(5);
        }
    }
    return "";
}

std::string library_key(const std::vector<std::string>& fields, const std::map<std::string, std::string>& libraries) {
    const auto read_group = record_read_group(fields);
    if (read_group.empty()) {
        return "unknown";
    }
    auto found = libraries.find(read_group);
    return found == libraries.end() ? read_group : found->second;
}

std::optional<DuplicateTemplate> make_single_duplicate_template(const std::vector<SamRecordLine>& records,
                                                               std::size_t record_index,
                                                               const std::map<std::string, std::string>& libraries) {
    const auto endpoint = duplicate_endpoint(records[record_index].fields);
    if (!endpoint) {
        return std::nullopt;
    }
    DuplicateTemplate duplicate_template;
    duplicate_template.record_indices.push_back(record_index);
    duplicate_template.key = library_key(records[record_index].fields, libraries) + "|SE|" + endpoint_key(*endpoint);
    duplicate_template.base_quality_score = base_quality_score(records[record_index].fields);
    duplicate_template.mapq_score = mapq_score(records[record_index].fields);
    duplicate_template.input_order = record_index;
    return duplicate_template;
}

std::optional<DuplicateTemplate> make_pair_duplicate_template(const std::vector<SamRecordLine>& records,
                                                             std::size_t first_index,
                                                             std::size_t second_index,
                                                             const std::map<std::string, std::string>& libraries) {
    auto first_endpoint = duplicate_endpoint(records[first_index].fields);
    auto second_endpoint = duplicate_endpoint(records[second_index].fields);
    if (!first_endpoint || !second_endpoint) {
        return std::nullopt;
    }
    auto first_key = endpoint_key(*first_endpoint);
    auto second_key = endpoint_key(*second_endpoint);
    if (second_key < first_key) {
        std::swap(first_key, second_key);
    }
    DuplicateTemplate duplicate_template;
    duplicate_template.record_indices.push_back(first_index);
    duplicate_template.record_indices.push_back(second_index);
    duplicate_template.key = library_key(records[first_index].fields, libraries) + "|PE|" + first_key + '|' + second_key;
    duplicate_template.base_quality_score = base_quality_score(records[first_index].fields) + base_quality_score(records[second_index].fields);
    duplicate_template.mapq_score = mapq_score(records[first_index].fields) + mapq_score(records[second_index].fields);
    duplicate_template.input_order = std::min(first_index, second_index);
    return duplicate_template;
}

std::string mark_duplicates_metrics_json(const MarkDuplicatesMetrics& metrics) {
    std::ostringstream out;
    out << "{\n";
    out << "  \"schema_version\": 1,\n";
    out << "  \"records\": " << metrics.records << ",\n";
    out << "  \"eligible_records\": " << metrics.eligible_records << ",\n";
    out << "  \"duplicate_sets\": " << metrics.duplicate_sets << ",\n";
    out << "  \"duplicates_marked\": " << metrics.duplicates_marked << ",\n";
    out << "  \"library_duplicates\": " << metrics.library_duplicates << ",\n";
    out << "  \"optical_duplicates\": " << metrics.optical_duplicates << ",\n";
    out << "  \"removed_duplicates\": " << metrics.removed_duplicates << "\n";
    out << "}\n";
    return out.str();
}

} // namespace

Result<SamValidationSummary> validate_sam_text(const std::filesystem::path& input,
                                                const ReferenceDictionary* reference_dict,
                                                bool summary_only) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    SamValidationSummary summary;
    std::ostringstream report;
    std::map<std::string, std::int64_t> header_lengths;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (line.empty()) {
            continue;
        }
        if (line[0] == '@') {
            if (starts_with(line, "@SQ")) {
                auto fields = split(line, '\t');
                std::string name;
                std::int64_t length = -1;
                for (const auto& field : fields) {
                    if (starts_with(field, "SN:")) {
                        name = field.substr(3);
                    } else if (starts_with(field, "LN:")) {
                        parse_i64(std::string_view(field).substr(3), length);
                    }
                }
                if (!name.empty() && length >= 0) {
                    header_lengths[name] = length;
                }
            }
            continue;
        }
        ++summary.records;
        auto fields = split(line, '\t');
        auto add_error = [&](const std::string& message) {
            ++summary.errors;
            if (!summary_only) {
                report << "line " << line_number << ": " << message << "\n";
            }
        };
        if (fields.size() < 11) {
            add_error("record has fewer than 11 SAM fields");
            continue;
        }
        std::int64_t pos = 0;
        if (!parse_i64(fields[3], pos) || pos < 0) {
            add_error("invalid POS");
        }
        if (!is_unmapped_flag(fields[1]) && fields[2] == "*") {
            add_error("mapped record has RNAME '*'");
        }
        auto cigar = parse_cigar(fields[5]);
        if (!cigar) {
            add_error(cigar.error().message);
            continue;
        }
        auto stats = cigar_stats(cigar.value());
        if (fields[9] != "*" && stats.query_length != static_cast<std::int64_t>(fields[9].size())) {
            add_error("CIGAR query length does not match SEQ length");
        }
        if (fields[2] != "*" && reference_dict != nullptr && reference_dict->lengths.find(fields[2]) == reference_dict->lengths.end()) {
            add_error("RNAME missing from reference dictionary: " + fields[2]);
        }
        auto header_found = header_lengths.find(fields[2]);
        const std::int64_t ref_length = header_found != header_lengths.end()
                                            ? header_found->second
                                            : (reference_dict != nullptr && reference_dict->lengths.count(fields[2]) != 0 ? reference_dict->lengths.at(fields[2]) : -1);
        if (pos > 0 && ref_length >= 0 && pos - 1 + stats.reference_span > ref_length) {
            add_error("alignment reference span exceeds contig length");
        }
    }
    report << "records=" << summary.records << " errors=" << summary.errors << "\n";
    summary.report = report.str();
    return summary;
}

Result<void> sort_sam_text(const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::string& order,
                           bool force) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::vector<std::string> headers;
    std::vector<SamRecordLine> records;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '@') {
            headers.push_back(line);
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        SamRecordLine record;
        record.line = line;
        record.fields = split(line, '\t');
        if (record.fields.size() >= 4) {
            parse_i64(record.fields[3], record.pos);
        }
        records.push_back(record);
    }
    if (order == "queryname") {
        std::stable_sort(records.begin(), records.end(), [](const SamRecordLine& a, const SamRecordLine& b) {
            return a.fields.empty() || b.fields.empty() ? a.line < b.line : a.fields[0] < b.fields[0];
        });
    } else if (order == "coordinate") {
        std::stable_sort(records.begin(), records.end(), [](const SamRecordLine& a, const SamRecordLine& b) {
            const auto ar = a.fields.size() > 2 ? a.fields[2] : "*";
            const auto br = b.fields.size() > 2 ? b.fields[2] : "*";
            if (ar != br) {
                return ar < br;
            }
            if (a.pos != b.pos) {
                return a.pos < b.pos;
            }
            return a.line < b.line;
        });
    } else {
        return make_error("--order must be coordinate or queryname");
    }
    std::ostringstream out;
    for (const auto& header : headers) {
        out << header << "\n";
    }
    for (const auto& record : records) {
        out << record.line << "\n";
    }
    return write_text_file_atomic(output, out.str(), force);
}

Result<void> clean_sam_text(const std::filesystem::path& input,
                            const std::filesystem::path& output,
                            const ReferenceDictionary* reference_dict,
                            bool force) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::ostringstream out;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (!line.empty() && line[0] == '@') {
            out << line << "\n";
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 11) {
            return make_error("record has fewer than 11 SAM fields on line " + std::to_string(line_number));
        }
        const auto flag = parse_flag_or_zero(fields[1]);
        if ((flag & 0x4) != 0) {
            fields[2] = "*";
            fields[3] = "0";
            fields[4] = "0";
            fields[5] = "*";
            fields[6] = "*";
            fields[7] = "0";
            fields[8] = "0";
            out << join_tab(fields) << "\n";
            continue;
        }
        std::int64_t pos = 0;
        std::string reason;
        auto cigar = parse_cigar(fields[5]);
        if (fields[2] == "*") {
            reason = "mapped_rname_missing";
        } else if (!parse_i64(fields[3], pos) || pos <= 0) {
            reason = "invalid_pos";
        } else if (!cigar) {
            reason = "invalid_cigar";
        } else if (fields[9] != "*" && cigar_stats(cigar.value()).query_length != static_cast<std::int64_t>(fields[9].size())) {
            reason = "query_length_mismatch";
        } else if (reference_dict != nullptr) {
            auto found = reference_dict->lengths.find(fields[2]);
            if (found == reference_dict->lengths.end()) {
                reason = "reference_conflict";
            } else if (pos - 1 + cigar_stats(cigar.value()).reference_span > found->second) {
                reason = "reference_overflow";
            }
        }
        if (!reason.empty()) {
            mark_record_unmapped(fields, reason);
        }
        out << join_tab(fields) << "\n";
    }
    return write_text_file_atomic(output, out.str(), force);
}

Result<void> reorder_sam_text(const std::filesystem::path& input,
                              const std::filesystem::path& output,
                              const std::filesystem::path& reference_fai,
                              bool allow_missing_contigs,
                              bool force) {
    auto contigs = read_ordered_fai(reference_fai);
    if (!contigs) {
        return contigs.error();
    }
    std::map<std::string, std::int64_t> contig_lengths;
    for (const auto& contig : contigs.value()) {
        contig_lengths[contig.name] = contig.length;
    }
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::vector<std::string> hd_headers;
    std::vector<std::string> other_headers;
    std::vector<std::string> records;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (!line.empty() && line[0] == '@') {
            if (starts_with(line, "@HD")) {
                hd_headers.push_back(line);
            } else if (!starts_with(line, "@SQ")) {
                other_headers.push_back(line);
            }
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 11) {
            return make_error("record has fewer than 11 SAM fields on line " + std::to_string(line_number));
        }
        if (fields[2] != "*" && contig_lengths.find(fields[2]) == contig_lengths.end() && !allow_missing_contigs) {
            return make_error("record references contig absent from reference FAI on line " + std::to_string(line_number) + ": " + fields[2]);
        }
        records.push_back(line);
    }
    std::ostringstream out;
    if (hd_headers.empty()) {
        out << "@HD\tVN:1.6\tSO:unknown\n";
    } else {
        for (const auto& header : hd_headers) {
            out << header << "\n";
        }
    }
    for (const auto& contig : contigs.value()) {
        out << "@SQ\tSN:" << contig.name << "\tLN:" << contig.length << "\n";
    }
    for (const auto& header : other_headers) {
        out << header << "\n";
    }
    for (const auto& record : records) {
        out << record << "\n";
    }
    return write_text_file_atomic(output, out.str(), force);
}

Result<void> fix_mate_sam_text(const std::filesystem::path& input,
                               const std::filesystem::path& output,
                               bool force) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::vector<std::string> headers;
    std::vector<SamRecordLine> records;
    std::map<std::string, std::vector<std::size_t>> by_qname;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (!line.empty() && line[0] == '@') {
            headers.push_back(line);
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 11) {
            return make_error("record has fewer than 11 SAM fields on line " + std::to_string(line_number));
        }
        SamRecordLine record;
        record.fields = std::move(fields);
        record.line = line;
        records.push_back(std::move(record));
        by_qname[records.back().fields[0]].push_back(records.size() - 1);
    }
    for (const auto& [qname, indices] : by_qname) {
        (void)qname;
        if (indices.size() != 2) {
            continue;
        }
        auto& first = records[indices[0]].fields;
        auto& second = records[indices[1]].fields;
        auto first_flag = parse_flag_or_zero(first[1]);
        auto second_flag = parse_flag_or_zero(second[1]);
        const bool first_mapped = record_is_mapped(first);
        const bool second_mapped = record_is_mapped(second);
        set_or_clear_flag(first_flag, 0x1, true);
        set_or_clear_flag(second_flag, 0x1, true);
        set_or_clear_flag(first_flag, 0x8, !second_mapped);
        set_or_clear_flag(second_flag, 0x8, !first_mapped);
        set_or_clear_flag(first_flag, 0x20, second_mapped && ((second_flag & 0x10) != 0));
        set_or_clear_flag(second_flag, 0x20, first_mapped && ((first_flag & 0x10) != 0));
        if ((first_flag & 0xC0) == 0 && (second_flag & 0xC0) == 0) {
            first_flag |= 0x40;
            second_flag |= 0x80;
        }
        first[1] = std::to_string(first_flag);
        second[1] = std::to_string(second_flag);
        first[6] = second_mapped ? (first[2] == second[2] ? "=" : second[2]) : "*";
        second[6] = first_mapped ? (first[2] == second[2] ? "=" : first[2]) : "*";
        first[7] = second_mapped ? second[3] : "0";
        second[7] = first_mapped ? first[3] : "0";
        first[8] = "0";
        second[8] = "0";
        if (first_mapped && second_mapped && first[2] == second[2]) {
            auto first_end = record_reference_end(first);
            auto second_end = record_reference_end(second);
            std::int64_t first_pos = 0;
            std::int64_t second_pos = 0;
            if (first_end && second_end && parse_i64(first[3], first_pos) && parse_i64(second[3], second_pos)) {
                const auto template_start = std::min(first_pos, second_pos);
                const auto template_end = std::max(*first_end, *second_end);
                const auto template_length = template_end - template_start + 1;
                const bool first_is_left = first_pos < second_pos || (first_pos == second_pos && indices[0] < indices[1]);
                first[8] = std::to_string(first_is_left ? template_length : -template_length);
                second[8] = std::to_string(first_is_left ? -template_length : template_length);
            }
        }
        add_mate_tags(first, second, second_mapped);
        add_mate_tags(second, first, first_mapped);
    }
    std::ostringstream out;
    for (const auto& header : headers) {
        out << header << "\n";
    }
    for (const auto& record : records) {
        out << join_tab(record.fields) << "\n";
    }
    return write_text_file_atomic(output, out.str(), force);
}

Result<MarkDuplicatesMetrics> mark_duplicates_sam_text(const std::filesystem::path& input,
                                                       const std::filesystem::path& output,
                                                       const std::filesystem::path* metrics_output,
                                                       bool remove_duplicates,
                                                       bool duplicate_type_tags,
                                                       bool force) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::vector<std::string> headers;
    std::vector<SamRecordLine> records;
    std::map<std::string, std::vector<std::size_t>> by_qname;
    MarkDuplicatesMetrics metrics;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '@') {
            headers.push_back(line);
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        SamRecordLine record;
        record.line = line;
        record.fields = split(line, '\t');
        records.push_back(std::move(record));
        ++metrics.records;
        if (!records.back().fields.empty()) {
            by_qname[records.back().fields[0]].push_back(records.size() - 1);
        }
    }

    const auto libraries = read_group_libraries(headers);
    std::vector<DuplicateTemplate> duplicate_templates;
    for (const auto& [qname, indices] : by_qname) {
        (void)qname;
        std::vector<std::size_t> eligible_indices;
        for (const auto index : indices) {
            if (duplicate_record_eligible(records[index].fields)) {
                eligible_indices.push_back(index);
                ++metrics.eligible_records;
            }
        }
        if (eligible_indices.size() == 2 &&
            (parse_flag_or_zero(records[eligible_indices[0]].fields[1]) & 0x1) != 0 &&
            (parse_flag_or_zero(records[eligible_indices[1]].fields[1]) & 0x1) != 0) {
            auto duplicate_template = make_pair_duplicate_template(records, eligible_indices[0], eligible_indices[1], libraries);
            if (duplicate_template) {
                duplicate_templates.push_back(*duplicate_template);
            }
            continue;
        }
        for (const auto index : eligible_indices) {
            auto duplicate_template = make_single_duplicate_template(records, index, libraries);
            if (duplicate_template) {
                duplicate_templates.push_back(*duplicate_template);
            }
        }
    }

    std::map<std::string, std::vector<std::size_t>> template_groups;
    for (std::size_t i = 0; i < duplicate_templates.size(); ++i) {
        template_groups[duplicate_templates[i].key].push_back(i);
    }

    std::vector<bool> duplicate_record(records.size(), false);
    for (const auto& [key, template_indices] : template_groups) {
        (void)key;
        if (template_indices.size() < 2) {
            continue;
        }
        ++metrics.duplicate_sets;
        auto best = template_indices.front();
        for (const auto candidate : template_indices) {
            const auto& candidate_template = duplicate_templates[candidate];
            const auto& best_template = duplicate_templates[best];
            if (candidate_template.base_quality_score > best_template.base_quality_score ||
                (candidate_template.base_quality_score == best_template.base_quality_score && candidate_template.mapq_score > best_template.mapq_score) ||
                (candidate_template.base_quality_score == best_template.base_quality_score && candidate_template.mapq_score == best_template.mapq_score &&
                 candidate_template.input_order < best_template.input_order)) {
                best = candidate;
            }
        }
        for (const auto duplicate_template_index : template_indices) {
            if (duplicate_template_index == best) {
                continue;
            }
            for (const auto record_index : duplicate_templates[duplicate_template_index].record_indices) {
                duplicate_record[record_index] = true;
            }
        }
    }

    std::ostringstream out;
    for (const auto& header : headers) {
        out << header << "\n";
    }
    for (std::size_t i = 0; i < records.size(); ++i) {
        auto& fields = records[i].fields;
        if (fields.size() >= 11 && duplicate_record_eligible(fields)) {
            auto flag = parse_flag_or_zero(fields[1]);
            const bool is_duplicate = duplicate_record[i];
            set_or_clear_flag(flag, 0x400, is_duplicate);
            fields[1] = std::to_string(flag);
            remove_tags_with_prefixes(fields, {"DT:Z:"});
            if (is_duplicate) {
                ++metrics.duplicates_marked;
                ++metrics.library_duplicates;
                if (duplicate_type_tags) {
                    fields.push_back("DT:Z:LB");
                }
                if (remove_duplicates) {
                    ++metrics.removed_duplicates;
                    continue;
                }
            }
            out << join_tab(fields) << "\n";
        } else {
            out << records[i].line << "\n";
        }
    }

    auto wrote = write_text_file_atomic(output, out.str(), force);
    if (!wrote) {
        return wrote.error();
    }
    if (metrics_output != nullptr) {
        auto wrote_metrics = write_text_file_atomic(*metrics_output, mark_duplicates_metrics_json(metrics), force);
        if (!wrote_metrics) {
            return wrote_metrics.error();
        }
    }
    return metrics;
}

Result<void> lift_sam_text(const ChainIndex& index,
                           const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::filesystem::path* unmapped_output,
                           const std::filesystem::path* metrics_output,
                           bool dry_run,
                           bool force,
                           bool reason_tags) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::ostringstream lifted_out;
    std::ostringstream rejected_out;
    int total = 0;
    int lifted = 0;
    int rejected = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '@') {
            if (!dry_run) {
                lifted_out << line << "\n";
                if (unmapped_output != nullptr) {
                    rejected_out << line << "\n";
                }
            }
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        ++total;
        auto fields = split(line, '\t');
        if (fields.size() < 11 || is_unmapped_flag(fields[1]) || fields[2] == "*") {
            ++rejected;
            if (!dry_run && unmapped_output != nullptr) {
                rejected_out << line << "\n";
            }
            continue;
        }
        std::int64_t pos = 0;
        if (!parse_i64(fields[3], pos) || pos <= 0) {
            ++rejected;
            if (!dry_run && unmapped_output != nullptr) {
                rejected_out << line << "\n";
            }
            continue;
        }
        auto mapped = map_position(index, fields[2], pos - 1);
        if (!mapped) {
            ++rejected;
            if (!dry_run && unmapped_output != nullptr) {
                rejected_out << line << "\n";
            }
            continue;
        }
        ++lifted;
        if (!dry_run) {
            fields[2] = mapped->contig;
            fields[3] = std::to_string(mapped->position + 1);
            if (reason_tags) {
                fields.push_back("BR:Z:lifted");
                fields.push_back(std::string("BC:Z:") + mapped->chain_id);
            }
            lifted_out << join_tab(fields) << "\n";
        }
    }
    if (!dry_run) {
        auto wrote = write_text_file_atomic(output, lifted_out.str(), force);
        if (!wrote) {
            return wrote;
        }
        if (unmapped_output != nullptr) {
            auto wrote_rejected = write_text_file_atomic(*unmapped_output, rejected_out.str(), force);
            if (!wrote_rejected) {
                return wrote_rejected;
            }
        }
    }
    if (metrics_output != nullptr) {
        auto wrote_metrics = write_text_file_atomic(*metrics_output, metrics_json(total, lifted, rejected), force);
        if (!wrote_metrics) {
            return wrote_metrics;
        }
    }
    return Result<void>();
}

Result<std::string> explain_sam_read(const ChainIndex& index,
                                     const std::filesystem::path& input,
                                     const std::string& qname) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '@') {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.empty() || fields[0] != qname) {
            continue;
        }
        std::ostringstream out;
        out << "qname: " << qname << "\n";
        out << "line: " << line_number << "\n";
        if (fields.size() < 11) {
            out << "reason: malformed_sam_record\n";
            return out.str();
        }
        out << "source: " << fields[2] << ':' << fields[3] << "\n";
        if (is_unmapped_flag(fields[1]) || fields[2] == "*") {
            out << "reason: pre_existing_unmapped\n";
            return out.str();
        }
        std::int64_t pos = 0;
        if (!parse_i64(fields[3], pos) || pos <= 0) {
            out << "reason: invalid_pos\n";
            return out.str();
        }
        auto mapped = map_position(index, fields[2], pos - 1);
        if (!mapped) {
            out << "reason: missing_mapping\n";
            return out.str();
        }
        out << "destination: " << mapped->contig << ':' << mapped->position + 1 << "\n";
        out << "strand: " << mapped->strand << "\n";
        out << "chain_id: " << mapped->chain_id << "\n";
        out << "reason: lifted\n";
        return out.str();
    }
    return make_error("QNAME not found: " + qname);
}

Result<void> create_sequence_dictionary_from_fai(const std::filesystem::path& fai_path,
                                                 const std::filesystem::path& output,
                                                 const std::string& assembly,
                                                 const std::string& species,
                                                 const std::string& uri,
                                                 bool force) {
    std::ifstream in(fai_path);
    if (!in) {
        return make_error("failed to open FAI: " + fai_path.string());
    }
    std::ostringstream out;
    out << "@HD\tVN:1.6\tSO:unsorted\n";
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 2) {
            return make_error("invalid FAI line " + std::to_string(line_number));
        }
        std::int64_t length = 0;
        if (!parse_i64(fields[1], length) || length < 0) {
            return make_error("invalid FAI length on line " + std::to_string(line_number));
        }
        out << "@SQ\tSN:" << fields[0] << "\tLN:" << length;
        if (!assembly.empty()) {
            out << "\tAS:" << assembly;
        }
        if (!species.empty()) {
            out << "\tSP:" << species;
        }
        if (!uri.empty()) {
            out << "\tUR:" << uri;
        }
        out << "\n";
    }
    return write_text_file_atomic(output, out.str(), force);
}

Result<void> replace_read_groups_sam_text(const std::filesystem::path& input,
                                          const std::filesystem::path& output,
                                          const ReadGroupOptions& read_group,
                                          bool force) {
    if (read_group.id.empty() || read_group.sample.empty()) {
        return make_error("read group --id and --sample are required");
    }
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    std::vector<std::string> headers;
    std::vector<std::string> records;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '@') {
            if (!starts_with(line, "@RG")) {
                headers.push_back(line);
            }
            continue;
        }
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() >= 11) {
            remove_read_group_tag(fields);
            fields.push_back("RG:Z:" + read_group.id);
            records.push_back(join_tab(fields));
        } else {
            records.push_back(line);
        }
    }
    std::ostringstream out;
    bool inserted = false;
    for (const auto& header : headers) {
        out << header << "\n";
        if (!inserted && starts_with(header, "@HD")) {
            out << read_group_header(read_group) << "\n";
            inserted = true;
        }
    }
    if (!inserted) {
        out << read_group_header(read_group) << "\n";
    }
    for (const auto& record : records) {
        out << record << "\n";
    }
    return write_text_file_atomic(output, out.str(), force);
}

Result<AlignmentSummary> collect_alignment_summary_text(const std::filesystem::path& input) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open SAM: " + input.string());
    }
    AlignmentSummary summary;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '@') {
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 11) {
            return make_error("record has fewer than 11 SAM fields on line " + std::to_string(line_number));
        }
        const auto flag = parse_flag_or_zero(fields[1]);
        ++summary.total_records;
        if ((flag & 0x4) != 0) {
            ++summary.unmapped_records;
        } else {
            ++summary.mapped_records;
            std::int64_t mapq = 0;
            if (parse_i64(fields[4], mapq)) {
                summary.mapq_sum += mapq;
            }
        }
        if ((flag & 0x100) != 0) {
            ++summary.secondary_records;
        }
        if ((flag & 0x800) != 0) {
            ++summary.supplementary_records;
        }
        if ((flag & 0x400) != 0) {
            ++summary.duplicate_records;
        }
        if (fields[9] != "*") {
            summary.total_bases += static_cast<std::int64_t>(fields[9].size());
        }
        auto cigar = parse_cigar(fields[5]);
        if (!cigar) {
            return cigar.error();
        }
        summary.aligned_bases += cigar_stats(cigar.value()).query_length;
    }
    return summary;
}

std::string format_alignment_summary(const AlignmentSummary& summary, bool as_json) {
    const double pct_mapped = summary.total_records == 0 ? 0.0 : static_cast<double>(summary.mapped_records) / static_cast<double>(summary.total_records);
    const double mean_mapq = summary.mapped_records == 0 ? 0.0 : static_cast<double>(summary.mapq_sum) / static_cast<double>(summary.mapped_records);
    std::ostringstream out;
    if (as_json) {
        out << "{\n";
        out << "  \"schema_version\": 1,\n";
        out << "  \"total_records\": " << summary.total_records << ",\n";
        out << "  \"mapped_records\": " << summary.mapped_records << ",\n";
        out << "  \"unmapped_records\": " << summary.unmapped_records << ",\n";
        out << "  \"secondary_records\": " << summary.secondary_records << ",\n";
        out << "  \"supplementary_records\": " << summary.supplementary_records << ",\n";
        out << "  \"duplicate_records\": " << summary.duplicate_records << ",\n";
        out << "  \"total_bases\": " << summary.total_bases << ",\n";
        out << "  \"aligned_bases\": " << summary.aligned_bases << ",\n";
        out << "  \"pct_mapped\": " << pct_mapped << ",\n";
        out << "  \"mean_mapq\": " << mean_mapq << "\n";
        out << "}\n";
        return out.str();
    }
    out << "TOTAL_RECORDS\t" << summary.total_records << "\n";
    out << "MAPPED_RECORDS\t" << summary.mapped_records << "\n";
    out << "UNMAPPED_RECORDS\t" << summary.unmapped_records << "\n";
    out << "SECONDARY_RECORDS\t" << summary.secondary_records << "\n";
    out << "SUPPLEMENTARY_RECORDS\t" << summary.supplementary_records << "\n";
    out << "DUPLICATE_RECORDS\t" << summary.duplicate_records << "\n";
    out << "TOTAL_BASES\t" << summary.total_bases << "\n";
    out << "ALIGNED_BASES\t" << summary.aligned_bases << "\n";
    out << "PCT_MAPPED\t" << pct_mapped << "\n";
    out << "MEAN_MAPQ\t" << mean_mapq << "\n";
    return out.str();
}

} // namespace burnham
