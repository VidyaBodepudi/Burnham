#include "burnham/hts.hpp"

#include "burnham/text.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

#if BURNHAM_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/sam.h>
#endif

namespace burnham {
namespace {

std::string lowercase_extension(const std::filesystem::path& path) {
    auto ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return ext;
}

std::string alignment_format_name(const std::filesystem::path& path) {
    const auto ext = lowercase_extension(path);
    if (ext == ".bam") {
        return "BAM";
    }
    if (ext == ".cram") {
        return "CRAM";
    }
    if (ext == ".sam") {
        return "SAM";
    }
    return "alignment";
}

} // namespace

HtslibStatus htslib_status() {
#if BURNHAM_HAS_HTSLIB
    HtslibStatus status;
    status.enabled = true;
    status.version = hts_version();
    status.message = std::string("enabled version ") + status.version;
    return status;
#else
    return HtslibStatus{false, "", "disabled; configure with BURNHAM_ENABLE_HTSLIB=ON and vcpkg to enable BAM/CRAM probing"};
#endif
}

bool is_binary_alignment_path(const std::filesystem::path& path) {
    const auto ext = lowercase_extension(path);
    return ext == ".bam" || ext == ".cram";
}

Result<void> require_text_alignment_path(const std::filesystem::path& path, const std::string& command_name) {
    if (!is_binary_alignment_path(path)) {
        return Result<void>();
    }
    return make_error(command_name + " currently transforms text SAM only; detected " + alignment_format_name(path) +
                      " input. Use validate-sam for the current htslib binary IO probe, or finish the binary command port.");
}

Result<HtsAlignmentSummary> validate_hts_alignment_input(const std::filesystem::path& input) {
    if (!is_binary_alignment_path(input)) {
        return make_error("validate_hts_alignment_input expects BAM or CRAM input: " + input.string());
    }
#if BURNHAM_HAS_HTSLIB
    samFile* file = sam_open(input.string().c_str(), "r");
    if (file == nullptr) {
        return make_error("failed to open " + alignment_format_name(input) + " with htslib: " + input.string());
    }
    sam_hdr_t* header = sam_hdr_read(file);
    if (header == nullptr) {
        sam_close(file);
        return make_error("failed to read alignment header with htslib: " + input.string());
    }
    bam1_t* record = bam_init1();
    if (record == nullptr) {
        sam_hdr_destroy(header);
        sam_close(file);
        return make_error("failed to allocate htslib BAM record");
    }

    HtsAlignmentSummary summary;
    summary.format = alignment_format_name(input);
    while (true) {
        const int read_result = sam_read1(file, header, record);
        if (read_result >= 0) {
            ++summary.records;
            continue;
        }
        if (read_result == -1) {
            break;
        }
        bam_destroy1(record);
        sam_hdr_destroy(header);
        sam_close(file);
        return make_error("failed while reading alignment records with htslib: " + input.string());
    }

    bam_destroy1(record);
    sam_hdr_destroy(header);
    sam_close(file);
    return summary;
#else
    return make_error("unsupported_format: " + alignment_format_name(input) +
                      " input requires htslib; rebuild with BURNHAM_ENABLE_HTSLIB=ON and vcpkg");
#endif
}

std::string format_hts_alignment_summary(const HtsAlignmentSummary& summary, bool as_json) {
    std::ostringstream out;
    if (as_json) {
        out << "{\n";
        out << "  \"schema_version\": 1,\n";
        out << "  \"format\": \"" << json_escape(summary.format) << "\",\n";
        out << "  \"records\": " << summary.records << "\n";
        out << "}\n";
        return out.str();
    }
    out << "format: " << summary.format << "\n";
    out << "records: " << summary.records << "\n";
    return out.str();
}

} // namespace burnham