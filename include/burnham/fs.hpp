#pragma once

#include "burnham/result.hpp"

#include <filesystem>
#include <string>

namespace burnham {

Result<std::string> read_text_file(const std::filesystem::path& path);
Result<void> write_text_file_atomic(const std::filesystem::path& path, const std::string& data, bool force);
Result<void> ensure_can_write_output(const std::filesystem::path& path, bool force);

} // namespace burnham
