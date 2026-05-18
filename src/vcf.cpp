#include "burnham/vcf.hpp"

#include "burnham/fs.hpp"
#include "burnham/text.hpp"

#include <fstream>
#include <sstream>

namespace burnham {

Result<void> lift_vcf_text(const ChainIndex& index,
                           const std::filesystem::path& input,
                           const std::filesystem::path& output,
                           const std::filesystem::path* rejected_output,
                           const std::filesystem::path* metrics_output,
                           bool dry_run,
                           bool force) {
    std::ifstream in(input);
    if (!in) {
        return make_error("failed to open VCF: " + input.string());
    }
    std::ostringstream lifted_out;
    std::ostringstream rejected_out;
    int total = 0;
    int lifted = 0;
    int rejected = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (starts_with(line, "#")) {
            if (!dry_run) {
                lifted_out << line << "\n";
                if (rejected_output != nullptr) {
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
        if (fields.size() < 8) {
            ++rejected;
            if (!dry_run && rejected_output != nullptr) {
                rejected_out << line << "\tBURNHAM_REJECT=malformed\n";
            }
            continue;
        }
        std::int64_t pos = 0;
        if (!parse_i64(fields[1], pos) || pos <= 0) {
            ++rejected;
            if (!dry_run && rejected_output != nullptr) {
                rejected_out << line << "\tBURNHAM_REJECT=invalid_pos\n";
            }
            continue;
        }
        auto mapped = map_position(index, fields[0], pos - 1);
        if (!mapped) {
            ++rejected;
            if (!dry_run && rejected_output != nullptr) {
                rejected_out << line << "\tBURNHAM_REJECT=missing_mapping\n";
            }
            continue;
        }
        ++lifted;
        if (!dry_run) {
            fields[0] = mapped->contig;
            fields[1] = std::to_string(mapped->position + 1);
            lifted_out << join_tab(fields) << "\n";
        }
    }
    if (!dry_run) {
        auto wrote = write_text_file_atomic(output, lifted_out.str(), force);
        if (!wrote) {
            return wrote;
        }
        if (rejected_output != nullptr) {
            auto wrote_rejected = write_text_file_atomic(*rejected_output, rejected_out.str(), force);
            if (!wrote_rejected) {
                return wrote_rejected;
            }
        }
    }
    if (metrics_output != nullptr) {
        std::ostringstream metrics;
        metrics << "{\n";
        metrics << "  \"schema_version\": 1,\n";
        metrics << "  \"records\": " << total << ",\n";
        metrics << "  \"lifted\": " << lifted << ",\n";
        metrics << "  \"rejected\": " << rejected << "\n";
        metrics << "}\n";
        auto wrote_metrics = write_text_file_atomic(*metrics_output, metrics.str(), force);
        if (!wrote_metrics) {
            return wrote_metrics;
        }
    }
    return Result<void>();
}

} // namespace burnham
