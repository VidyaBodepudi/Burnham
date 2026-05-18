#include "burnham/chain_index.hpp"

#include "burnham/fs.hpp"
#include "burnham/text.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace burnham {
namespace {

struct ChainHeader {
    std::string t_name;
    std::int64_t t_size = 0;
    char t_strand = '+';
    std::int64_t t_start = 0;
    std::int64_t t_end = 0;
    std::string q_name;
    std::int64_t q_size = 0;
    char q_strand = '+';
    std::int64_t q_start = 0;
    std::int64_t q_end = 0;
    std::string id;
};

Result<ChainHeader> parse_header(const std::string& line) {
    auto fields = split_whitespace(line);
    if (fields.size() != 13 || fields[0] != "chain") {
        return make_error("invalid chain header: " + line);
    }
    ChainHeader header;
    header.t_name = fields[2];
    header.t_strand = fields[4].empty() ? '+' : fields[4][0];
    header.q_name = fields[7];
    header.q_strand = fields[9].empty() ? '+' : fields[9][0];
    header.id = fields[12];
    if (!parse_i64(fields[3], header.t_size) || !parse_i64(fields[5], header.t_start) ||
        !parse_i64(fields[6], header.t_end) || !parse_i64(fields[8], header.q_size) ||
        !parse_i64(fields[10], header.q_start) || !parse_i64(fields[11], header.q_end)) {
        return make_error("invalid numeric value in chain header: " + line);
    }
    if (header.t_strand != '+' || (header.q_strand != '+' && header.q_strand != '-')) {
        return make_error("unsupported chain strand in header: " + line);
    }
    return header;
}

bool dict_has(const ReferenceDictionary* dict, const std::string& contig, std::int64_t expected_size) {
    if (dict == nullptr) {
        return true;
    }
    auto found = dict->lengths.find(contig);
    return found != dict->lengths.end() && found->second == expected_size;
}

Result<void> validate_header_dicts(const ChainHeader& header,
                                   SourceSide source_side,
                                   const ReferenceDictionary* source_dict,
                                   const ReferenceDictionary* dest_dict) {
    const auto& source_name = source_side == SourceSide::Query ? header.q_name : header.t_name;
    const auto source_size = source_side == SourceSide::Query ? header.q_size : header.t_size;
    const auto& dest_name = source_side == SourceSide::Query ? header.t_name : header.q_name;
    const auto dest_size = source_side == SourceSide::Query ? header.t_size : header.q_size;
    if (!dict_has(source_dict, source_name, source_size)) {
        return make_error("source dictionary mismatch for contig: " + source_name);
    }
    if (!dict_has(dest_dict, dest_name, dest_size)) {
        return make_error("destination dictionary mismatch for contig: " + dest_name);
    }
    return Result<void>();
}

ChainBlock make_block(const ChainHeader& header, SourceSide source_side, std::int64_t t_pos, std::int64_t q_pos, std::int64_t size) {
    const std::int64_t t_start = t_pos;
    const std::int64_t t_end = t_pos + size;
    std::int64_t q_start = q_pos;
    std::int64_t q_end = q_pos + size;
    if (header.q_strand == '-') {
        q_start = header.q_size - (q_pos + size);
        q_end = header.q_size - q_pos;
    }

    ChainBlock block;
    block.chain_id = header.id;
    block.strand = header.q_strand == '-' ? '-' : '+';
    if (source_side == SourceSide::Query) {
        block.source_contig = header.q_name;
        block.source_start = q_start;
        block.source_end = q_end;
        block.dest_contig = header.t_name;
        block.dest_start = t_start;
        block.dest_end = t_end;
    } else {
        block.source_contig = header.t_name;
        block.source_start = t_start;
        block.source_end = t_end;
        block.dest_contig = header.q_name;
        block.dest_start = q_start;
        block.dest_end = q_end;
    }
    return block;
}

} // namespace

std::string source_side_name(SourceSide side) {
    return side == SourceSide::Query ? "query" : "target";
}

Result<SourceSide> parse_source_side(const std::string& value) {
    if (value == "query") {
        return SourceSide::Query;
    }
    if (value == "target") {
        return SourceSide::Target;
    }
    return make_error("--source-side must be query or target");
}

Result<ReferenceDictionary> read_fai(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        return make_error("failed to open FAI: " + path.string());
    }
    ReferenceDictionary dict;
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
        dict.lengths[fields[0]] = length;
    }
    return dict;
}

Result<ChainIndex> build_chain_index(const std::filesystem::path& chain_path,
                                     SourceSide source_side,
                                     const ReferenceDictionary* source_dict,
                                     const ReferenceDictionary* dest_dict) {
    std::ifstream in(chain_path);
    if (!in) {
        return make_error("failed to open chain file: " + chain_path.string());
    }

    ChainIndex index;
    index.source_side = source_side;
    std::string line;
    std::optional<ChainHeader> current;
    std::int64_t t_pos = 0;
    std::int64_t q_pos = 0;
    int line_number = 0;

    while (std::getline(in, line)) {
        ++line_number;
        const auto trimmed = trim(line);
        if (trimmed.empty()) {
            current.reset();
            continue;
        }
        if (starts_with(trimmed, "chain ")) {
            auto header = parse_header(trimmed);
            if (!header) {
                return header.error();
            }
            auto dict_check = validate_header_dicts(header.value(), source_side, source_dict, dest_dict);
            if (!dict_check) {
                return dict_check.error();
            }
            current = header.value();
            t_pos = current->t_start;
            q_pos = current->q_start;
            continue;
        }
        if (!current) {
            return make_error("chain block data before header on line " + std::to_string(line_number));
        }
        auto fields = split_whitespace(trimmed);
        if (fields.size() != 1 && fields.size() != 3) {
            return make_error("invalid chain block line " + std::to_string(line_number));
        }
        std::int64_t size = 0;
        if (!parse_i64(fields[0], size) || size <= 0) {
            return make_error("invalid chain block size on line " + std::to_string(line_number));
        }
        index.blocks.push_back(make_block(*current, source_side, t_pos, q_pos, size));
        t_pos += size;
        q_pos += size;
        if (fields.size() == 3) {
            std::int64_t dt = 0;
            std::int64_t dq = 0;
            if (!parse_i64(fields[1], dt) || !parse_i64(fields[2], dq) || dt < 0 || dq < 0) {
                return make_error("invalid chain gap on line " + std::to_string(line_number));
            }
            t_pos += dt;
            q_pos += dq;
        } else {
            current.reset();
        }
    }
    std::sort(index.blocks.begin(), index.blocks.end(), [](const ChainBlock& a, const ChainBlock& b) {
        if (a.source_contig != b.source_contig) {
            return a.source_contig < b.source_contig;
        }
        return a.source_start < b.source_start;
    });
    return index;
}

Result<void> write_chain_index(const ChainIndex& index, const std::filesystem::path& path, bool force) {
    std::ostringstream out;
    out << "BURNHAM_CHAIN_INDEX\t1\n";
    out << "source_side\t" << source_side_name(index.source_side) << "\n";
    for (const auto& block : index.blocks) {
        out << "block\t" << block.source_contig << '\t' << block.source_start << '\t' << block.source_end << '\t'
            << block.dest_contig << '\t' << block.dest_start << '\t' << block.dest_end << '\t'
            << block.strand << '\t' << block.chain_id << "\n";
    }
    return write_text_file_atomic(path, out.str(), force);
}

Result<ChainIndex> read_chain_index(const std::filesystem::path& path) {
    auto data = read_text_file(path);
    if (!data) {
        return data.error();
    }
    std::istringstream in(data.value());
    ChainIndex index;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (line_number == 1) {
            if (fields.size() != 2 || fields[0] != "BURNHAM_CHAIN_INDEX" || fields[1] != "1") {
                return make_error("unsupported chain index format");
            }
            continue;
        }
        if (fields[0] == "source_side" && fields.size() == 2) {
            auto side = parse_source_side(fields[1]);
            if (!side) {
                return side.error();
            }
            index.source_side = side.value();
            continue;
        }
        if (fields[0] == "block" && fields.size() == 9) {
            ChainBlock block;
            block.source_contig = fields[1];
            block.dest_contig = fields[4];
            block.strand = fields[7].empty() ? '+' : fields[7][0];
            block.chain_id = fields[8];
            if (!parse_i64(fields[2], block.source_start) || !parse_i64(fields[3], block.source_end) ||
                !parse_i64(fields[5], block.dest_start) || !parse_i64(fields[6], block.dest_end)) {
                return make_error("invalid numeric value in chain index line " + std::to_string(line_number));
            }
            index.blocks.push_back(block);
            continue;
        }
        return make_error("invalid chain index line " + std::to_string(line_number));
    }
    return index;
}

std::optional<MappedPosition> map_position(const ChainIndex& index, const std::string& contig, std::int64_t position) {
    for (const auto& block : index.blocks) {
        if (block.source_contig != contig || position < block.source_start || position >= block.source_end) {
            continue;
        }
        const auto offset = position - block.source_start;
        MappedPosition mapped;
        mapped.contig = block.dest_contig;
        mapped.strand = block.strand;
        mapped.chain_id = block.chain_id;
        mapped.position = block.strand == '+' ? block.dest_start + offset : block.dest_end - 1 - offset;
        return mapped;
    }
    return std::nullopt;
}

std::optional<MappedInterval> map_interval(const ChainIndex& index, const std::string& contig, std::int64_t start, std::int64_t end) {
    if (start < 0 || end <= start) {
        return std::nullopt;
    }
    for (const auto& block : index.blocks) {
        if (block.source_contig != contig || start < block.source_start || end > block.source_end) {
            continue;
        }
        const auto start_offset = start - block.source_start;
        const auto last_offset = (end - 1) - block.source_start;
        const auto mapped_start = block.strand == '+' ? block.dest_start + start_offset : block.dest_end - 1 - start_offset;
        const auto mapped_last = block.strand == '+' ? block.dest_start + last_offset : block.dest_end - 1 - last_offset;
        MappedInterval mapped;
        mapped.contig = block.dest_contig;
        mapped.strand = block.strand;
        mapped.chain_id = block.chain_id;
        mapped.start = std::min(mapped_start, mapped_last);
        mapped.end = std::max(mapped_start, mapped_last) + 1;
        return mapped;
    }
    return std::nullopt;
}

std::string inspect_chain_index_text(const ChainIndex& index, bool as_json) {
    std::map<std::string, std::int64_t> coverage;
    int plus = 0;
    int minus = 0;
    std::set<std::string> dest_contigs;
    for (const auto& block : index.blocks) {
        coverage[block.source_contig] += block.source_end - block.source_start;
        dest_contigs.insert(block.dest_contig);
        if (block.strand == '-') {
            ++minus;
        } else {
            ++plus;
        }
    }
    std::ostringstream out;
    if (as_json) {
        out << "{\n";
        out << "  \"format\": \"burnham-chain-index\",\n";
        out << "  \"version\": 1,\n";
        out << "  \"source_side\": \"" << source_side_name(index.source_side) << "\",\n";
        out << "  \"blocks\": " << index.blocks.size() << ",\n";
        out << "  \"plus_blocks\": " << plus << ",\n";
        out << "  \"minus_blocks\": " << minus << ",\n";
        out << "  \"source_contigs\": " << coverage.size() << ",\n";
        out << "  \"dest_contigs\": " << dest_contigs.size() << "\n";
        out << "}\n";
        return out.str();
    }
    out << "Burnham chain index v1\n";
    out << "source_side: " << source_side_name(index.source_side) << "\n";
    out << "blocks: " << index.blocks.size() << "\n";
    out << "strand_blocks: +=" << plus << " -= " << minus << "\n";
    out << "source_contigs: " << coverage.size() << "\n";
    for (const auto& [contig, bases] : coverage) {
        out << "  " << contig << " covered_bases=" << bases << "\n";
    }
    out << "dest_contigs: " << dest_contigs.size() << "\n";
    return out.str();
}

} // namespace burnham
