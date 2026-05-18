#pragma once

#include "burnham/chain_index.hpp"
#include "burnham/result.hpp"

#include <filesystem>

namespace burnham {

Result<void> lift_vcf_text(const ChainIndex& index,
                           const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::filesystem::path* rejected_output,
                           const std::filesystem::path* metrics_output,
                           bool dry_run,
                           bool force);

} // namespace burnham
