#include "burnham/text.hpp"

#include <charconv>
#include <cctype>
#include <sstream>

namespace burnham {

std::string trim(std::string_view input) {
    std::size_t first = 0;
    while (first < input.size() && std::isspace(static_cast<unsigned char>(input[first])) != 0) {
        ++first;
    }
    std::size_t last = input.size();
    while (last > first && std::isspace(static_cast<unsigned char>(input[last - 1])) != 0) {
        --last;
    }
    return std::string(input.substr(first, last - first));
}

std::vector<std::string> split(std::string_view input, char delimiter) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= input.size()) {
        const std::size_t pos = input.find(delimiter, start);
        if (pos == std::string_view::npos) {
            fields.emplace_back(input.substr(start));
            break;
        }
        fields.emplace_back(input.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

std::vector<std::string> split_whitespace(std::string_view input) {
    std::vector<std::string> fields;
    std::istringstream stream{std::string(input)};
    std::string field;
    while (stream >> field) {
        fields.push_back(field);
    }
    return fields;
}

bool starts_with(std::string_view input, std::string_view prefix) {
    return input.size() >= prefix.size() && input.substr(0, prefix.size()) == prefix;
}

bool parse_i64(std::string_view input, std::int64_t& value) {
    if (input.empty()) {
        return false;
    }
    const char* begin = input.data();
    const char* end = input.data() + input.size();
    auto [ptr, ec] = std::from_chars(begin, end, value);
    return ec == std::errc() && ptr == end;
}

std::string join_tab(const std::vector<std::string>& fields) {
    std::ostringstream out;
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (i != 0) {
            out << '\t';
        }
        out << fields[i];
    }
    return out.str();
}

std::string json_escape(std::string_view input) {
    std::ostringstream out;
    for (char ch : input) {
        switch (ch) {
        case '\\': out << "\\\\"; break;
        case '"': out << "\\\""; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << ch; break;
        }
    }
    return out.str();
}

} // namespace burnham
