#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace burnham {

std::string trim(std::string_view input);
std::vector<std::string> split(std::string_view input, char delimiter);
std::vector<std::string> split_whitespace(std::string_view input);
bool starts_with(std::string_view input, std::string_view prefix);
bool parse_i64(std::string_view input, std::int64_t& value);
std::string join_tab(const std::vector<std::string>& fields);
std::string json_escape(std::string_view input);

} // namespace burnham
