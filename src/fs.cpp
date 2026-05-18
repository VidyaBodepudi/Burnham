#include "burnham/fs.hpp"

#include <fstream>
#include <sstream>

namespace burnham {

Result<std::string> read_text_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return make_error("failed to open input file: " + path.string());
    }
    std::ostringstream data;
    data << in.rdbuf();
    return data.str();
}

Result<void> ensure_can_write_output(const std::filesystem::path& path, bool force) {
    if (path.empty()) {
        return make_error("output path is empty");
    }
    if (std::filesystem::exists(path) && !force) {
        return make_error("output already exists, use --force: " + path.string());
    }
    const auto parent = path.parent_path();
    if (!parent.empty() && !std::filesystem::exists(parent)) {
        return make_error("output directory does not exist: " + parent.string());
    }
    return Result<void>();
}

Result<void> write_text_file_atomic(const std::filesystem::path& path, const std::string& data, bool force) {
    auto can_write = ensure_can_write_output(path, force);
    if (!can_write) {
        return can_write;
    }

    const auto temp = path.string() + ".tmp";
    {
        std::ofstream out(temp, std::ios::binary | std::ios::trunc);
        if (!out) {
            return make_error("failed to create temporary output: " + temp);
        }
        out << data;
        if (!out) {
            return make_error("failed while writing temporary output: " + temp);
        }
    }

    std::error_code error;
    if (force && std::filesystem::exists(path)) {
        std::filesystem::remove(path, error);
        if (error) {
            std::filesystem::remove(temp);
            return make_error("failed to replace existing output: " + path.string());
        }
    }
    std::filesystem::rename(temp, path, error);
    if (error) {
        std::filesystem::remove(temp);
        return make_error("failed to publish output: " + path.string());
    }
    return Result<void>();
}

} // namespace burnham
