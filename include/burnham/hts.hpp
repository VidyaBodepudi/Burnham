#pragma once

#include "burnham/result.hpp"

#include <cstdint>
#include <filesystem>
#include <string>

#ifndef BURNHAM_HAS_HTSLIB
#define BURNHAM_HAS_HTSLIB 0
#endif

namespace burnham {

struct HtslibStatus {
    bool enabled = false;
    std::string version;
    std::string message;
};

struct HtsAlignmentSummary {
    std::int64_t records = 0;
    std::string format;
};

HtslibStatus htslib_status();
bool is_binary_alignment_path(const std::filesystem::path& path);
Result<void> require_text_alignment_path(const std::filesystem::path& path, const std::string& command_name);
Result<HtsAlignmentSummary> validate_hts_alignment_input(const std::filesystem::path& input);
std::string format_hts_alignment_summary(const HtsAlignmentSummary& summary, bool as_json);

} // namespace burnham