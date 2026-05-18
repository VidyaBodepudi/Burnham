#pragma once

#include "burnham/result.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace burnham {

enum class SampleHaplotype { All, First, Second };

struct SampleVariant {
    std::string contig;
    std::int64_t position = 0;
    std::string id;
    std::string ref;
    std::string alt;
    std::int64_t delta = 0;
};

struct SampleIndex {
    std::string sample = "all";
    SampleHaplotype haplotype = SampleHaplotype::All;
    std::vector<SampleVariant> variants;
};

struct SampleMappingResult {
    bool mapped = false;
    std::string contig;
    std::int64_t position = 0;
    std::string reason;
};

Result<SampleHaplotype> parse_sample_haplotype(const std::string& value);
std::string sample_haplotype_name(SampleHaplotype haplotype);
Result<SampleIndex> build_sample_index(const std::filesystem::path& vcf_path,
                                       const std::string& sample,
                                       SampleHaplotype haplotype);
Result<void> write_sample_index(const SampleIndex& index, const std::filesystem::path& path, bool force);
Result<SampleIndex> read_sample_index(const std::filesystem::path& path);
SampleMappingResult map_sample_position(const SampleIndex& index, const std::string& contig, std::int64_t position);
std::optional<SampleMappingResult> map_sample_position_optional(const SampleIndex& index, const std::string& contig, std::int64_t position);
std::string inspect_sample_index_text(const SampleIndex& index, bool as_json);
Result<void> lift_sam_with_sample_index_text(const SampleIndex& index,
                                             const std::filesystem::path& input,
                                             const std::filesystem::path& output,
                                             const std::filesystem::path* unmapped_output,
                                             const std::filesystem::path* metrics_output,
                                             bool dry_run,
                                             bool force,
                                             bool reason_tags);

} // namespace burnham
