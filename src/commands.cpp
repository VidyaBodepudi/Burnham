#include "burnham/commands.hpp"

#include "burnham/chain_index.hpp"
#include "burnham/fs.hpp"
#include "burnham/hts.hpp"
#include "burnham/sample_index.hpp"
#include "burnham/sam.hpp"
#include "burnham/text.hpp"
#include "burnham/vcf.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace burnham {
namespace {

struct ParsedArgs {
    std::map<std::string, std::string> values;
    std::map<std::string, bool> flags;
};

ParsedArgs parse_args(int argc, char** argv, int start) {
    ParsedArgs parsed;
    for (int i = start; i < argc; ++i) {
        std::string arg = argv[i];
        if (!starts_with(arg, "--")) {
            parsed.flags[arg] = true;
            continue;
        }
        if (arg == "--force" || arg == "--dry-run" || arg == "--json" || arg == "--summary" || arg == "--reason-tags" ||
            arg == "--allow-missing-contigs" || arg == "--remove-duplicates" || arg == "--duplicate-type-tags") {
            parsed.flags[arg] = true;
            continue;
        }
        if (i + 1 >= argc) {
            parsed.values[arg] = "";
            continue;
        }
        parsed.values[arg] = argv[++i];
    }
    return parsed;
}

std::optional<std::string> get_value(const ParsedArgs& args, const std::string& key) {
    auto found = args.values.find(key);
    if (found == args.values.end() || found->second.empty()) {
        return std::nullopt;
    }
    return found->second;
}

bool has_flag(const ParsedArgs& args, const std::string& key) {
    auto found = args.flags.find(key);
    return found != args.flags.end() && found->second;
}

Result<std::filesystem::path> require_path(const ParsedArgs& args, const std::string& key) {
    auto value = get_value(args, key);
    if (!value) {
        return make_error("missing required option: " + key);
    }
    return std::filesystem::path(*value);
}

void print_help(std::ostream& out) {
    out << "Burnham 0.1.0\n";
    out << "Commands:\n";
    out << "  doctor\n";
    out << "  index-chain --chain input.chain --output output.bci --source-side query|target\n";
    out << "  index-sample --vcf input.vcf --output output.bsi [--sample SAMPLE] [--haplotype all|0|1]\n";
    out << "  inspect-index --index output.bci [--json]\n";
    out << "  inspect-sample-index --index output.bsi [--json]\n";
    out << "  lift-bed --index output.bci --input input.bed --output lifted.bed\n";
    out << "  validate-sam --input input.sam|input.bam|input.cram [--reference-fai ref.fai]\n";
    out << "  sort --input input.sam --output sorted.sam --order coordinate|queryname\n";
    out << "  dict --reference-fai ref.fai --output ref.dict\n";
    out << "  read-groups --input input.sam --output output.sam --id ID --sample SAMPLE\n";
    out << "  alignment-summary --input input.sam [--output metrics.txt] [--json]\n";
    out << "  clean --input input.sam --output clean.sam [--reference-fai ref.fai]\n";
    out << "  reorder --input input.sam --output reordered.sam --reference-fai ref.fai\n";
    out << "  fix-mate --input input.sam --output fixed.sam\n";
    out << "  mark-dup --input input.sam --output marked.sam [--metrics-file metrics.json]\n";
    out << "  lift-sam --index output.bci --input input.sam --output lifted.sam\n";
    out << "  lift-sam --sample-index output.bsi --input input.sam --output lifted.sam\n";
    out << "  explain-read --index output.bci --input input.sam --qname read-name\n";
    out << "  lift-vcf-chain --index output.bci --input input.vcf --output lifted.vcf\n";
}

Result<ChainIndex> load_index_from_args(const ParsedArgs& args) {
    auto index_path = require_path(args, "--index");
    if (!index_path) {
        return index_path.error();
    }
    return read_chain_index(index_path.value());
}

Result<void> command_index_chain(const ParsedArgs& args) {
    auto chain_path = require_path(args, "--chain");
    if (!chain_path) {
        return chain_path.error();
    }
    auto output_path = require_path(args, "--output");
    if (!output_path) {
        return output_path.error();
    }
    auto side_text = get_value(args, "--source-side");
    if (!side_text) {
        return make_error("missing required option: --source-side query|target");
    }
    auto side = parse_source_side(*side_text);
    if (!side) {
        return side.error();
    }
    std::optional<ReferenceDictionary> source_dict;
    std::optional<ReferenceDictionary> dest_dict;
    if (auto source_fai = get_value(args, "--source-fai")) {
        auto dict = read_fai(*source_fai);
        if (!dict) {
            return dict.error();
        }
        source_dict = dict.value();
    }
    if (auto dest_fai = get_value(args, "--dest-fai")) {
        auto dict = read_fai(*dest_fai);
        if (!dict) {
            return dict.error();
        }
        dest_dict = dict.value();
    }
    auto index = build_chain_index(chain_path.value(), side.value(), source_dict ? &*source_dict : nullptr, dest_dict ? &*dest_dict : nullptr);
    if (!index) {
        return index.error();
    }
    return write_chain_index(index.value(), output_path.value(), has_flag(args, "--force"));
}

Result<void> command_lift_bed(const ParsedArgs& args) {
    auto index = load_index_from_args(args);
    if (!index) {
        return index.error();
    }
    auto input_path = require_path(args, "--input");
    if (!input_path) {
        return input_path.error();
    }
    auto output_path = require_path(args, "--output");
    if (!output_path) {
        return output_path.error();
    }
    std::ifstream in(input_path.value());
    if (!in) {
        return make_error("failed to open BED: " + input_path.value().string());
    }
    std::ostringstream lifted_out;
    std::ostringstream unmapped_out;
    int total = 0;
    int lifted = 0;
    int unmapped = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (trim(line).empty() || starts_with(line, "#") || starts_with(line, "track") || starts_with(line, "browser")) {
            if (!has_flag(args, "--dry-run")) {
                lifted_out << line << "\n";
            }
            continue;
        }
        ++total;
        auto fields = split(line, '\t');
        if (fields.size() < 3) {
            ++unmapped;
            if (!has_flag(args, "--dry-run")) {
                unmapped_out << line << "\tinvalid_interval\n";
            }
            continue;
        }
        std::int64_t start = 0;
        std::int64_t end = 0;
        if (!parse_i64(fields[1], start) || !parse_i64(fields[2], end) || end <= start) {
            ++unmapped;
            if (!has_flag(args, "--dry-run")) {
                unmapped_out << line << "\tinvalid_interval\n";
            }
            continue;
        }
        auto mapped = map_interval(index.value(), fields[0], start, end);
        if (!mapped) {
            ++unmapped;
            if (!has_flag(args, "--dry-run")) {
                unmapped_out << line << "\tmissing_mapping\n";
            }
            continue;
        }
        ++lifted;
        if (!has_flag(args, "--dry-run")) {
            fields[0] = mapped->contig;
            fields[1] = std::to_string(mapped->start);
            fields[2] = std::to_string(mapped->end);
            lifted_out << join_tab(fields) << "\n";
        }
    }
    const bool force = has_flag(args, "--force");
    if (!has_flag(args, "--dry-run")) {
        auto wrote = write_text_file_atomic(output_path.value(), lifted_out.str(), force);
        if (!wrote) {
            return wrote;
        }
        if (auto unmapped_path = get_value(args, "--unmapped-output")) {
            auto wrote_unmapped = write_text_file_atomic(*unmapped_path, unmapped_out.str(), force);
            if (!wrote_unmapped) {
                return wrote_unmapped;
            }
        }
    }
    if (auto metrics_path = get_value(args, "--metrics-file")) {
        std::ostringstream metrics;
        metrics << "{\n";
        metrics << "  \"schema_version\": 1,\n";
        metrics << "  \"records\": " << total << ",\n";
        metrics << "  \"lifted\": " << lifted << ",\n";
        metrics << "  \"unmapped\": " << unmapped << "\n";
        metrics << "}\n";
        auto wrote_metrics = write_text_file_atomic(*metrics_path, metrics.str(), force);
        if (!wrote_metrics) {
            return wrote_metrics;
        }
    }
    return Result<void>();
}

Result<void> command_index_sample(const ParsedArgs& args) {
    auto vcf_path = require_path(args, "--vcf");
    if (!vcf_path) {
        return vcf_path.error();
    }
    auto output_path = require_path(args, "--output");
    if (!output_path) {
        return output_path.error();
    }
    const auto sample = get_value(args, "--sample").value_or("all");
    auto haplotype = parse_sample_haplotype(get_value(args, "--haplotype").value_or("all"));
    if (!haplotype) {
        return haplotype.error();
    }
    auto index = build_sample_index(vcf_path.value(), sample, haplotype.value());
    if (!index) {
        return index.error();
    }
    return write_sample_index(index.value(), output_path.value(), has_flag(args, "--force"));
}

} // namespace

int run_cli(int argc, char** argv, std::ostream& out, std::ostream& err) {
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        print_help(out);
        return argc < 2 ? 1 : 0;
    }
    const std::string command = argv[1];
    const auto args = parse_args(argc, argv, 2);
    auto fail = [&](const Error& error) {
        err << "error: " << error.message << "\n";
        return 1;
    };

    if (command == "doctor") {
        const auto hts_status = htslib_status();
        out << "Burnham doctor\n";
        out << "cxx_core: ok\n";
        out << "htslib: " << hts_status.message << "\n";
        if (auto index_path = get_value(args, "--chain-index")) {
            auto index = read_chain_index(*index_path);
            if (!index) {
                return fail(index.error());
            }
            out << "chain_index: ok blocks=" << index.value().blocks.size() << "\n";
        }
        if (auto sample_index_path = get_value(args, "--sample-index")) {
            auto index = read_sample_index(*sample_index_path);
            if (!index) {
                return fail(index.error());
            }
            out << "sample_index: ok variants=" << index.value().variants.size() << " sample=" << index.value().sample << "\n";
        }
        return 0;
    }
    if (command == "index-chain") {
        auto result = command_index_chain(args);
        return result ? 0 : fail(result.error());
    }
    if (command == "index-sample") {
        auto result = command_index_sample(args);
        return result ? 0 : fail(result.error());
    }
    if (command == "inspect-index") {
        auto index = load_index_from_args(args);
        if (!index) {
            return fail(index.error());
        }
        out << inspect_chain_index_text(index.value(), has_flag(args, "--json"));
        return 0;
    }
    if (command == "inspect-sample-index") {
        auto index_path = require_path(args, "--index");
        if (!index_path) {
            return fail(index_path.error());
        }
        auto index = read_sample_index(index_path.value());
        if (!index) {
            return fail(index.error());
        }
        out << inspect_sample_index_text(index.value(), has_flag(args, "--json"));
        return 0;
    }
    if (command == "lift-bed") {
        auto result = command_lift_bed(args);
        return result ? 0 : fail(result.error());
    }
    if (command == "validate-sam") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        if (is_binary_alignment_path(input_path.value())) {
            auto result = validate_hts_alignment_input(input_path.value());
            if (!result) {
                return fail(result.error());
            }
            out << format_hts_alignment_summary(result.value(), has_flag(args, "--json"));
            return 0;
        }
        std::optional<ReferenceDictionary> reference_dict;
        if (auto reference_fai = get_value(args, "--reference-fai")) {
            auto dict = read_fai(*reference_fai);
            if (!dict) {
                return fail(dict.error());
            }
            reference_dict = dict.value();
        }
        auto result = validate_sam_text(input_path.value(), reference_dict ? &*reference_dict : nullptr, has_flag(args, "--summary"));
        if (!result) {
            return fail(result.error());
        }
        out << result.value().report;
        return result.value().errors == 0 ? 0 : 2;
    }
    if (command == "sort") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        const auto order = get_value(args, "--order").value_or("coordinate");
        auto result = sort_sam_text(input_path.value(), output_path.value(), order, has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "dict") {
        auto fai_path = require_path(args, "--reference-fai");
        if (!fai_path) {
            return fail(fai_path.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        auto result = create_sequence_dictionary_from_fai(fai_path.value(), output_path.value(), get_value(args, "--assembly").value_or(""),
                                                         get_value(args, "--species").value_or(""), get_value(args, "--uri").value_or(""),
                                                         has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "read-groups") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        ReadGroupOptions read_group;
        read_group.id = get_value(args, "--id").value_or("");
        read_group.sample = get_value(args, "--sample").value_or("");
        read_group.library = get_value(args, "--library").value_or("");
        read_group.platform = get_value(args, "--platform").value_or("");
        read_group.platform_unit = get_value(args, "--platform-unit").value_or("");
        auto result = replace_read_groups_sam_text(input_path.value(), output_path.value(), read_group, has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "alignment-summary") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto summary = collect_alignment_summary_text(input_path.value());
        if (!summary) {
            return fail(summary.error());
        }
        auto formatted = format_alignment_summary(summary.value(), has_flag(args, "--json"));
        if (auto output_path = get_value(args, "--output")) {
            auto wrote = write_text_file_atomic(*output_path, formatted, has_flag(args, "--force"));
            return wrote ? 0 : fail(wrote.error());
        }
        out << formatted;
        return 0;
    }
    if (command == "clean") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        std::optional<ReferenceDictionary> reference_dict;
        if (auto reference_fai = get_value(args, "--reference-fai")) {
            auto dict = read_fai(*reference_fai);
            if (!dict) {
                return fail(dict.error());
            }
            reference_dict = dict.value();
        }
        auto result = clean_sam_text(input_path.value(), output_path.value(), reference_dict ? &*reference_dict : nullptr, has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "reorder") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        auto fai_path = require_path(args, "--reference-fai");
        if (!fai_path) {
            return fail(fai_path.error());
        }
        auto result = reorder_sam_text(input_path.value(), output_path.value(), fai_path.value(), has_flag(args, "--allow-missing-contigs"),
                                       has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "fix-mate") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        auto result = fix_mate_sam_text(input_path.value(), output_path.value(), has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "mark-dup") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        auto metrics_text = get_value(args, "--metrics-file");
        std::optional<std::filesystem::path> metrics_path;
        if (metrics_text) {
            metrics_path = *metrics_text;
        }
        auto result = mark_duplicates_sam_text(input_path.value(), output_path.value(), metrics_path ? &*metrics_path : nullptr,
                                               has_flag(args, "--remove-duplicates"), has_flag(args, "--duplicate-type-tags"),
                                               has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }
    if (command == "lift-sam") {
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        auto unmapped_text = get_value(args, "--unmapped-output");
        auto metrics_text = get_value(args, "--metrics-file");
        std::optional<std::filesystem::path> unmapped_path;
        std::optional<std::filesystem::path> metrics_path;
        if (unmapped_text) {
            unmapped_path = *unmapped_text;
        }
        if (metrics_text) {
            metrics_path = *metrics_text;
        }
        const auto sample_index_text = get_value(args, "--sample-index");
        const auto chain_index_text = get_value(args, "--index");
        if (sample_index_text && chain_index_text) {
            return fail(make_error("provide either --index or --sample-index, not both"));
        }
        if (sample_index_text) {
            auto sample_index = read_sample_index(*sample_index_text);
            if (!sample_index) {
                return fail(sample_index.error());
            }
            auto result = lift_sam_with_sample_index_text(sample_index.value(), input_path.value(), output_path.value(), unmapped_path ? &*unmapped_path : nullptr,
                                                          metrics_path ? &*metrics_path : nullptr, has_flag(args, "--dry-run"), has_flag(args, "--force"),
                                                          has_flag(args, "--reason-tags"));
            return result ? 0 : fail(result.error());
        }
        auto index = load_index_from_args(args);
        if (!index) {
            return fail(index.error());
        }
        auto result = lift_sam_text(index.value(), input_path.value(), output_path.value(), unmapped_path ? &*unmapped_path : nullptr,
                                    metrics_path ? &*metrics_path : nullptr, has_flag(args, "--dry-run"), has_flag(args, "--force"),
                                    has_flag(args, "--reason-tags"));
        return result ? 0 : fail(result.error());
    }
    if (command == "explain-read") {
        auto index = load_index_from_args(args);
        if (!index) {
            return fail(index.error());
        }
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto text_guard = require_text_alignment_path(input_path.value(), command);
        if (!text_guard) {
            return fail(text_guard.error());
        }
        auto qname = get_value(args, "--qname");
        if (!qname) {
            return fail(make_error("missing required option: --qname"));
        }
        auto result = explain_sam_read(index.value(), input_path.value(), *qname);
        if (!result) {
            return fail(result.error());
        }
        out << result.value();
        return 0;
    }
    if (command == "lift-vcf-chain") {
        auto index = load_index_from_args(args);
        if (!index) {
            return fail(index.error());
        }
        auto input_path = require_path(args, "--input");
        if (!input_path) {
            return fail(input_path.error());
        }
        auto output_path = require_path(args, "--output");
        if (!output_path) {
            return fail(output_path.error());
        }
        auto rejected_text = get_value(args, "--rejected-output");
        auto metrics_text = get_value(args, "--metrics-file");
        std::optional<std::filesystem::path> rejected_path;
        std::optional<std::filesystem::path> metrics_path;
        if (rejected_text) {
            rejected_path = *rejected_text;
        }
        if (metrics_text) {
            metrics_path = *metrics_text;
        }
        auto result = lift_vcf_text(index.value(), input_path.value(), output_path.value(), rejected_path ? &*rejected_path : nullptr,
                                    metrics_path ? &*metrics_path : nullptr, has_flag(args, "--dry-run"), has_flag(args, "--force"));
        return result ? 0 : fail(result.error());
    }

    err << "error: unknown command: " << command << "\n";
    print_help(err);
    return 1;
}

} // namespace burnham
