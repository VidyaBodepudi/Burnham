#pragma once

#include "burnham/chain_index.hpp"
#include "burnham/result.hpp"

#include <filesystem>
#include <cstdint>
#include <string>

namespace burnham {

struct SamValidationSummary {
    int records = 0;
    int errors = 0;
    std::string report;
};

struct ReadGroupOptions {
    std::string id;
    std::string sample;
    std::string library;
    std::string platform;
    std::string platform_unit;
};

struct AlignmentSummary {
    std::int64_t total_records = 0;
    std::int64_t mapped_records = 0;
    std::int64_t unmapped_records = 0;
    std::int64_t secondary_records = 0;
    std::int64_t supplementary_records = 0;
    std::int64_t duplicate_records = 0;
    std::int64_t total_bases = 0;
    std::int64_t aligned_bases = 0;
    std::int64_t mapq_sum = 0;
};

struct MarkDuplicatesMetrics {
    std::int64_t records = 0;
    std::int64_t eligible_records = 0;
    std::int64_t duplicate_sets = 0;
    std::int64_t duplicates_marked = 0;
    std::int64_t library_duplicates = 0;
    std::int64_t optical_duplicates = 0;
    std::int64_t removed_duplicates = 0;
};

Result<SamValidationSummary> validate_sam_text(const std::filesystem::path& input,
                                                const ReferenceDictionary* reference_dict,
                                                bool summary_only);
Result<void> sort_sam_text(const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::string& order,
                           bool force);
Result<void> clean_sam_text(const std::filesystem::path& input,
                            const std::filesystem::path& output,
                            const ReferenceDictionary* reference_dict,
                            bool force);
Result<void> reorder_sam_text(const std::filesystem::path& input,
                              const std::filesystem::path& output,
                              const std::filesystem::path& reference_fai,
                              bool allow_missing_contigs,
                              bool force);
Result<void> fix_mate_sam_text(const std::filesystem::path& input,
                               const std::filesystem::path& output,
                               bool force);
Result<MarkDuplicatesMetrics> mark_duplicates_sam_text(const std::filesystem::path& input,
                                                       const std::filesystem::path& output,
                                                       const std::filesystem::path* metrics_output,
                                                       bool remove_duplicates,
                                                       bool duplicate_type_tags,
                                                       bool force);
Result<void> lift_sam_text(const ChainIndex& index,
                           const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::filesystem::path* unmapped_output,
                           const std::filesystem::path* metrics_output,
                           bool dry_run,
                           bool force,
                           bool reason_tags);
Result<std::string> explain_sam_read(const ChainIndex& index,
                                     const std::filesystem::path& input,
                                     const std::string& qname);
Result<void> create_sequence_dictionary_from_fai(const std::filesystem::path& fai_path,
                                                 const std::filesystem::path& output,
                                                 const std::string& assembly,
                                                 const std::string& species,
                                                 const std::string& uri,
                                                 bool force);
Result<void> replace_read_groups_sam_text(const std::filesystem::path& input,
                                          const std::filesystem::path& output,
                                          const ReadGroupOptions& read_group,
                                          bool force);
Result<AlignmentSummary> collect_alignment_summary_text(const std::filesystem::path& input);
std::string format_alignment_summary(const AlignmentSummary& summary, bool as_json);

} // namespace burnham
