#include "burnham/sample_index.hpp"

#include "burnham/fs.hpp"
#include "burnham/text.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>

namespace burnham {
namespace {

bool has_symbolic_or_multi_alt(const std::string& alt) {
    return alt.empty() || alt.find(',') != std::string::npos || alt.find('<') != std::string::npos || alt.find('>') != std::string::npos;
}

std::optional<std::size_t> genotype_field_index(const std::string& format) {
    auto fields = split(format, ':');
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (fields[i] == "GT") {
            return i;
        }
    }
    return std::nullopt;
}

std::vector<std::string> genotype_alleles(const std::string& genotype) {
    std::string normalized = genotype;
    std::replace(normalized.begin(), normalized.end(), '|', '/');
    return split(normalized, '/');
}

bool genotype_selects_alt(const std::string& format, const std::string& sample_value, SampleHaplotype haplotype) {
    auto gt_index = genotype_field_index(format);
    if (!gt_index) {
        return false;
    }
    auto sample_fields = split(sample_value, ':');
    if (*gt_index >= sample_fields.size()) {
        return false;
    }
    auto alleles = genotype_alleles(sample_fields[*gt_index]);
    if (haplotype == SampleHaplotype::All) {
        return std::find(alleles.begin(), alleles.end(), "1") != alleles.end();
    }
    const std::size_t allele_index = haplotype == SampleHaplotype::First ? 0 : 1;
    return allele_index < alleles.size() && alleles[allele_index] == "1";
}

bool is_snv(const SampleVariant& variant) {
    return variant.ref.size() == 1 && variant.alt.size() == 1;
}

bool is_insertion(const SampleVariant& variant) {
    return variant.alt.size() > variant.ref.size();
}

bool is_deletion(const SampleVariant& variant) {
    return variant.ref.size() > variant.alt.size();
}

bool is_mnp(const SampleVariant& variant) {
    return variant.ref.size() == variant.alt.size() && variant.ref.size() > 1;
}

std::string sample_metrics_json(int total, int lifted, int rejected, int deleted_by_variant) {
    std::ostringstream out;
    out << "{\n";
    out << "  \"schema_version\": 1,\n";
    out << "  \"engine\": \"sample\",\n";
    out << "  \"records\": " << total << ",\n";
    out << "  \"sample_lifted\": " << lifted << ",\n";
    out << "  \"rejected\": " << rejected << ",\n";
    out << "  \"deleted_by_variant\": " << deleted_by_variant << "\n";
    out << "}\n";
    return out.str();
}

bool is_unmapped_sam_flag(const std::string& flag_text) {
    std::int64_t flag = 0;
    return parse_i64(flag_text, flag) && ((flag & 0x4) != 0);
}

} // namespace

Result<SampleHaplotype> parse_sample_haplotype(const std::string& value) {
    if (value == "all") {
        return SampleHaplotype::All;
    }
    if (value == "0" || value == "first") {
        return SampleHaplotype::First;
    }
    if (value == "1" || value == "second") {
        return SampleHaplotype::Second;
    }
    return make_error("--haplotype must be all, 0, first, 1, or second");
}

std::string sample_haplotype_name(SampleHaplotype haplotype) {
    switch (haplotype) {
    case SampleHaplotype::First: return "0";
    case SampleHaplotype::Second: return "1";
    case SampleHaplotype::All: return "all";
    }
    return "all";
}

Result<SampleIndex> build_sample_index(const std::filesystem::path& vcf_path,
                                       const std::string& sample,
                                       SampleHaplotype haplotype) {
    std::ifstream in(vcf_path);
    if (!in) {
        return make_error("failed to open VCF: " + vcf_path.string());
    }

    SampleIndex index;
    index.sample = sample.empty() ? "all" : sample;
    index.haplotype = haplotype;
    std::optional<std::size_t> sample_column;
    bool saw_chrom_header = false;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (trim(line).empty()) {
            continue;
        }
        if (starts_with(line, "##")) {
            continue;
        }
        if (starts_with(line, "#CHROM")) {
            saw_chrom_header = true;
            auto fields = split(line, '\t');
            if (index.sample != "all") {
                for (std::size_t i = 9; i < fields.size(); ++i) {
                    if (fields[i] == index.sample) {
                        sample_column = i;
                        break;
                    }
                }
                if (!sample_column) {
                    return make_error("sample not found in VCF header: " + index.sample);
                }
            }
            continue;
        }
        auto fields = split(line, '\t');
        if (fields.size() < 8) {
            return make_error("malformed VCF record on line " + std::to_string(line_number));
        }
        if (!saw_chrom_header) {
            return make_error("VCF #CHROM header is required before records");
        }
        if (has_symbolic_or_multi_alt(fields[4])) {
            continue;
        }
        if (sample_column) {
            if (fields.size() <= *sample_column || fields.size() < 10 || !genotype_selects_alt(fields[8], fields[*sample_column], haplotype)) {
                continue;
            }
        }
        std::int64_t pos = 0;
        if (!parse_i64(fields[1], pos) || pos <= 0) {
            return make_error("invalid VCF POS on line " + std::to_string(line_number));
        }
        SampleVariant variant;
        variant.contig = fields[0];
        variant.position = pos - 1;
        variant.id = fields[2] == "." ? "" : fields[2];
        variant.ref = fields[3];
        variant.alt = fields[4];
        variant.delta = static_cast<std::int64_t>(variant.alt.size()) - static_cast<std::int64_t>(variant.ref.size());
        if (variant.ref.empty() || variant.alt.empty()) {
            continue;
        }
        index.variants.push_back(std::move(variant));
    }
    std::sort(index.variants.begin(), index.variants.end(), [](const SampleVariant& left, const SampleVariant& right) {
        if (left.contig != right.contig) {
            return left.contig < right.contig;
        }
        return left.position < right.position;
    });
    return index;
}

Result<void> write_sample_index(const SampleIndex& index, const std::filesystem::path& path, bool force) {
    std::ostringstream out;
    out << "BURNHAM_SAMPLE_INDEX\t1\n";
    out << "sample\t" << index.sample << "\n";
    out << "haplotype\t" << sample_haplotype_name(index.haplotype) << "\n";
    for (const auto& variant : index.variants) {
        out << "variant\t" << variant.contig << '\t' << variant.position << '\t' << variant.ref << '\t'
            << variant.alt << '\t' << variant.delta << '\t' << variant.id << "\n";
    }
    return write_text_file_atomic(path, out.str(), force);
}

Result<SampleIndex> read_sample_index(const std::filesystem::path& path) {
    auto data = read_text_file(path);
    if (!data) {
        return data.error();
    }
    SampleIndex index;
    std::istringstream in(data.value());
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (trim(line).empty()) {
            continue;
        }
        auto fields = split(line, '\t');
        if (line_number == 1) {
            if (fields.size() != 2 || fields[0] != "BURNHAM_SAMPLE_INDEX" || fields[1] != "1") {
                return make_error("unsupported sample index format");
            }
            continue;
        }
        if (fields[0] == "sample" && fields.size() == 2) {
            index.sample = fields[1];
            continue;
        }
        if (fields[0] == "haplotype" && fields.size() == 2) {
            auto haplotype = parse_sample_haplotype(fields[1]);
            if (!haplotype) {
                return haplotype.error();
            }
            index.haplotype = haplotype.value();
            continue;
        }
        if (fields[0] == "variant" && fields.size() == 7) {
            SampleVariant variant;
            variant.contig = fields[1];
            variant.ref = fields[3];
            variant.alt = fields[4];
            variant.id = fields[6];
            if (!parse_i64(fields[2], variant.position) || !parse_i64(fields[5], variant.delta)) {
                return make_error("invalid sample index numeric value on line " + std::to_string(line_number));
            }
            index.variants.push_back(std::move(variant));
            continue;
        }
        return make_error("invalid sample index line " + std::to_string(line_number));
    }
    return index;
}

SampleMappingResult map_sample_position(const SampleIndex& index, const std::string& contig, std::int64_t position) {
    if (position < 0) {
        return SampleMappingResult{false, contig, position, "invalid_pos"};
    }
    std::int64_t shift = 0;
    for (const auto& variant : index.variants) {
        if (variant.contig < contig) {
            continue;
        }
        if (variant.contig > contig) {
            break;
        }
        if (position < variant.position) {
            break;
        }
        const auto ref_length = static_cast<std::int64_t>(variant.ref.size());
        const auto alt_length = static_cast<std::int64_t>(variant.alt.size());
        const auto variant_end = variant.position + ref_length;
        if (position > variant.position && position < variant_end && alt_length < ref_length) {
            return SampleMappingResult{false, contig, position, "deleted_by_variant"};
        }
        if (position >= variant_end) {
            shift += variant.delta;
        }
    }
    return SampleMappingResult{true, contig, position + shift, "sample_lifted"};
}

std::optional<SampleMappingResult> map_sample_position_optional(const SampleIndex& index, const std::string& contig, std::int64_t position) {
    auto mapped = map_sample_position(index, contig, position);
    if (!mapped.mapped) {
        return std::nullopt;
    }
    return mapped;
}

std::string inspect_sample_index_text(const SampleIndex& index, bool as_json) {
    int snvs = 0;
    int insertions = 0;
    int deletions = 0;
    int mnps = 0;
    std::map<std::string, int> contigs;
    for (const auto& variant : index.variants) {
        ++contigs[variant.contig];
        if (is_snv(variant)) {
            ++snvs;
        } else if (is_insertion(variant)) {
            ++insertions;
        } else if (is_deletion(variant)) {
            ++deletions;
        } else if (is_mnp(variant)) {
            ++mnps;
        }
    }
    std::ostringstream out;
    if (as_json) {
        out << "{\n";
        out << "  \"format\": \"burnham-sample-index\",\n";
        out << "  \"version\": 1,\n";
        out << "  \"sample\": \"" << json_escape(index.sample) << "\",\n";
        out << "  \"haplotype\": \"" << sample_haplotype_name(index.haplotype) << "\",\n";
        out << "  \"variants\": " << index.variants.size() << ",\n";
        out << "  \"contigs\": " << contigs.size() << ",\n";
        out << "  \"snvs\": " << snvs << ",\n";
        out << "  \"insertions\": " << insertions << ",\n";
        out << "  \"deletions\": " << deletions << ",\n";
        out << "  \"mnps\": " << mnps << "\n";
        out << "}\n";
        return out.str();
    }
    out << "Burnham sample index v1\n";
    out << "sample: " << index.sample << "\n";
    out << "haplotype: " << sample_haplotype_name(index.haplotype) << "\n";
    out << "variants: " << index.variants.size() << "\n";
    out << "contigs: " << contigs.size() << "\n";
    out << "snvs: " << snvs << "\n";
    out << "insertions: " << insertions << "\n";
    out << "deletions: " << deletions << "\n";
    out << "mnps: " << mnps << "\n";
    return out.str();
}

Result<void> lift_sam_with_sample_index_text(const SampleIndex& index,
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
    int deleted_by_variant = 0;
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
        if (fields.size() < 11 || is_unmapped_sam_flag(fields[1]) || fields[2] == "*") {
            ++rejected;
            if (!dry_run && unmapped_output != nullptr) {
                rejected_out << line << "\tBR:Z:pre_existing_unmapped\n";
            }
            continue;
        }
        std::int64_t pos = 0;
        if (!parse_i64(fields[3], pos) || pos <= 0) {
            ++rejected;
            if (!dry_run && unmapped_output != nullptr) {
                rejected_out << line << "\tBR:Z:invalid_pos\n";
            }
            continue;
        }
        auto mapped = map_sample_position(index, fields[2], pos - 1);
        if (!mapped.mapped) {
            ++rejected;
            if (mapped.reason == "deleted_by_variant") {
                ++deleted_by_variant;
            }
            if (!dry_run && unmapped_output != nullptr) {
                rejected_out << line << "\tBR:Z:" << mapped.reason << "\n";
            }
            continue;
        }
        ++lifted;
        if (!dry_run) {
            fields[2] = mapped.contig;
            fields[3] = std::to_string(mapped.position + 1);
            if (reason_tags) {
                fields.push_back("BR:Z:" + mapped.reason);
                fields.push_back("BS:Z:" + index.sample);
                fields.push_back("BH:Z:" + sample_haplotype_name(index.haplotype));
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
        auto wrote_metrics = write_text_file_atomic(*metrics_output, sample_metrics_json(total, lifted, rejected, deleted_by_variant), force);
        if (!wrote_metrics) {
            return wrote_metrics;
        }
    }
    return Result<void>();
}

} // namespace burnham
