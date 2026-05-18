#include "burnham/chain_index.hpp"
#include "burnham/cigar.hpp"
#include "burnham/fs.hpp"
#include "burnham/hts.hpp"
#include "burnham/sample_index.hpp"
#include "burnham/sam.hpp"
#include "burnham/vcf.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

using namespace burnham;

std::filesystem::path temp_root() {
    auto path = std::filesystem::temp_directory_path() / "burnham_tests";
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    return path;
}

void write_file(const std::filesystem::path& path, const std::string& data) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << data;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "test failed: " << message << "\n";
        std::exit(1);
    }
}

void test_cigar() {
    auto cigar = parse_cigar("5S10M2I3D4=");
    require(cigar.ok(), "parse valid cigar");
    auto stats = cigar_stats(cigar.value());
    require(stats.query_length == 21, "cigar query length");
    require(stats.reference_span == 17, "cigar reference span");
    require(!parse_cigar("0M").ok(), "reject zero cigar op");
}

void test_hts_status_and_binary_guards() {
    auto status = htslib_status();
    require(!status.message.empty(), "htslib status message");
    require(is_binary_alignment_path("sample.bam"), "detect bam path");
    require(is_binary_alignment_path("sample.cram"), "detect cram path");
    require(!is_binary_alignment_path("sample.sam"), "sam path remains text");

    auto guard = require_text_alignment_path("reads.bam", "sort");
    require(!guard.ok(), "binary transform guard rejects bam");
    require(guard.error().message.find("text SAM only") != std::string::npos, "binary transform guard message");

    auto probe = validate_hts_alignment_input("missing.bam");
    require(!probe.ok(), "missing bam probe fails cleanly");
    if (!status.enabled) {
        require(probe.error().message.find("requires htslib") != std::string::npos, "disabled htslib probe message");
    }
}

void test_chain_index_and_interval_mapping() {
    auto root = temp_root();
    auto chain = root / "simple.chain";
    write_file(chain,
               "chain 1 chrDst 1000 + 100 220 chrSrc 1000 + 10 130 7\n"
               "50 5 5\n"
               "65\n");
    auto index = build_chain_index(chain, SourceSide::Query, nullptr, nullptr);
    require(index.ok(), "build chain index");
    require(index.value().blocks.size() == 2, "chain block count");
    auto mapped = map_interval(index.value(), "chrSrc", 20, 30);
    require(mapped.has_value(), "map interval");
    require(mapped->contig == "chrDst", "mapped contig");
    require(mapped->start == 110 && mapped->end == 120, "mapped interval coordinates");
    require(!map_interval(index.value(), "chrSrc", 55, 70).has_value(), "split interval rejected");

    auto index_path = root / "simple.bci";
    auto wrote = write_chain_index(index.value(), index_path, true);
    require(wrote.ok(), "write index");
    auto loaded = read_chain_index(index_path);
    require(loaded.ok(), "read index");
    require(loaded.value().blocks.size() == 2, "loaded index block count");
}

void test_fai_validation() {
    auto root = temp_root();
    auto fai = root / "ref.fai";
    write_file(fai, "chrSrc\t1000\t0\t80\t81\n");
    auto dict = read_fai(fai);
    require(dict.ok(), "read fai");
    require(dict.value().lengths.at("chrSrc") == 1000, "fai contig length");
}

void test_sam_validate_sort_lift() {
    auto root = temp_root();
    auto chain = root / "simple.chain";
    write_file(chain,
               "chain 1 chrDst 1000 + 100 200 chrSrc 1000 + 0 100 1\n"
               "100\n");
    auto index = build_chain_index(chain, SourceSide::Query, nullptr, nullptr);
    require(index.ok(), "build sam chain index");

    auto sam = root / "input.sam";
    write_file(sam,
               "@HD\tVN:1.6\tSO:unknown\n"
               "@SQ\tSN:chrSrc\tLN:1000\n"
               "readB\t0\tchrSrc\t20\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"
               "readA\t0\tchrSrc\t10\t60\t4M\t*\t0\t0\tACGT\t!!!!\n");
    auto validation = validate_sam_text(sam, nullptr, false);
    require(validation.ok(), "validate sam text");
    require(validation.value().errors == 0, "sam validation errors");

    auto sorted = root / "sorted.sam";
    auto sorted_result = sort_sam_text(sam, sorted, "queryname", true);
    require(sorted_result.ok(), "sort sam text");
    auto sorted_data = read_text_file(sorted);
    require(sorted_data.ok(), "read sorted sam");
    require(sorted_data.value().find("readA") < sorted_data.value().find("readB"), "queryname sort order");

    auto lifted = root / "lifted.sam";
    auto metrics = root / "sam.metrics.json";
    auto lifted_result = lift_sam_text(index.value(), sam, lifted, nullptr, &metrics, false, true, true);
    require(lifted_result.ok(), "lift sam text");
    auto lifted_data = read_text_file(lifted);
    require(lifted_data.ok(), "read lifted sam");
    require(lifted_data.value().find("chrDst\t110") != std::string::npos, "lifted sam coordinate");
    require(lifted_data.value().find("BR:Z:lifted") != std::string::npos, "lifted sam reason tag");

    auto explanation = explain_sam_read(index.value(), sam, "readA");
    require(explanation.ok(), "explain read");
    require(explanation.value().find("destination: chrDst:110") != std::string::npos, "explain read destination");
}

void test_reverse_chain_mapping() {
    auto root = temp_root();
    auto chain = root / "reverse.chain";
    write_file(chain,
               "chain 1 chrDst 1000 + 100 110 chrSrc 1000 - 0 10 9\n"
               "10\n");
    auto index = build_chain_index(chain, SourceSide::Query, nullptr, nullptr);
    require(index.ok(), "build reverse chain index");
    auto mapped = map_interval(index.value(), "chrSrc", 990, 1000);
    require(mapped.has_value(), "map reverse interval");
    require(mapped->contig == "chrDst", "reverse mapped contig");
    require(mapped->start == 100 && mapped->end == 110, "reverse mapped coordinates");
    require(mapped->strand == '-', "reverse mapped strand");
}

void test_vcf_lift() {
    auto root = temp_root();
    auto chain = root / "simple.chain";
    write_file(chain,
               "chain 1 chrDst 1000 + 100 200 chrSrc 1000 + 0 100 1\n"
               "100\n");
    auto index = build_chain_index(chain, SourceSide::Query, nullptr, nullptr);
    require(index.ok(), "build vcf chain index");

    auto vcf = root / "input.vcf";
    write_file(vcf,
               "##fileformat=VCFv4.2\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "chrSrc\t11\t.\tA\tC\t.\tPASS\t.\n");
    auto lifted = root / "lifted.vcf";
    auto result = lift_vcf_text(index.value(), vcf, lifted, nullptr, nullptr, false, true);
    require(result.ok(), "lift vcf text");
    auto lifted_data = read_text_file(lifted);
    require(lifted_data.ok(), "read lifted vcf");
    require(lifted_data.value().find("chrDst\t111") != std::string::npos, "lifted vcf coordinate");
}

void test_sample_index_and_sam_lift() {
    auto root = temp_root();
    auto vcf = root / "sample.vcf";
    write_file(vcf,
               "##fileformat=VCFv4.2\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
               "chr1\t11\tins1\tA\tATG\t.\tPASS\t.\tGT\t0/1\t0/0\n"
               "chr1\t21\tdel1\tAT\tA\t.\tPASS\t.\tGT\t1/1\t0/0\n"
               "chr1\t31\tskip1\tC\tG\t.\tPASS\t.\tGT\t0/0\t1/1\n"
               "chr1\t41\tsymbolic\tN\t<DEL>\t.\tPASS\t.\tGT\t1/1\t1/1\n");
    auto index = build_sample_index(vcf, "S1", SampleHaplotype::All);
    require(index.ok(), "build sample index");
    require(index.value().variants.size() == 2, "sample variant count");

    auto shifted = map_sample_position(index.value(), "chr1", 14);
    require(shifted.mapped && shifted.position == 16, "sample insertion shift");
    auto deleted = map_sample_position(index.value(), "chr1", 21);
    require(!deleted.mapped && deleted.reason == "deleted_by_variant", "sample deletion rejection");
    auto after_deletion = map_sample_position(index.value(), "chr1", 24);
    require(after_deletion.mapped && after_deletion.position == 25, "sample deletion downstream shift");

    auto index_path = root / "sample.bsi";
    auto wrote = write_sample_index(index.value(), index_path, true);
    require(wrote.ok(), "write sample index");
    auto loaded = read_sample_index(index_path);
    require(loaded.ok(), "read sample index");
    require(loaded.value().sample == "S1", "sample index metadata");
    auto inspection = inspect_sample_index_text(loaded.value(), true);
    require(inspection.find("\"variants\": 2") != std::string::npos, "sample inspect variant count");

    auto sam = root / "sample.sam";
    write_file(sam,
               "@HD\tVN:1.6\tSO:unknown\n"
               "@SQ\tSN:chr1\tLN:1000\n"
               "readShift\t0\tchr1\t15\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"
               "readDeleted\t0\tchr1\t22\t60\t4M\t*\t0\t0\tACGT\t!!!!\n");
    auto lifted = root / "sample_lifted.sam";
    auto rejected = root / "sample_rejected.sam";
    auto metrics = root / "sample.metrics.json";
    auto lifted_result = lift_sam_with_sample_index_text(loaded.value(), sam, lifted, &rejected, &metrics, false, true, true);
    require(lifted_result.ok(), "sample lift sam text");
    auto lifted_data = read_text_file(lifted);
    require(lifted_data.ok(), "read sample lifted sam");
    require(lifted_data.value().find("readShift\t0\tchr1\t17") != std::string::npos, "sample lifted coordinate");
    require(lifted_data.value().find("BS:Z:S1") != std::string::npos, "sample lifted tag");
    auto rejected_data = read_text_file(rejected);
    require(rejected_data.ok(), "read sample rejected sam");
    require(rejected_data.value().find("BR:Z:deleted_by_variant") != std::string::npos, "sample rejected reason");
    auto metrics_data = read_text_file(metrics);
    require(metrics_data.ok(), "read sample metrics");
    require(metrics_data.value().find("\"sample_lifted\": 1") != std::string::npos, "sample metrics lifted count");
}

void test_phase3_sam_utilities() {
    auto root = temp_root();
    auto fai = root / "ref.fai";
    write_file(fai,
               "chr2\t200\t0\t80\t81\n"
               "chr1\t100\t210\t80\t81\n");
    auto dict = root / "ref.dict";
    auto dict_result = create_sequence_dictionary_from_fai(fai, dict, "asm1", "human", "file://ref.fa", true);
    require(dict_result.ok(), "create sequence dictionary");
    auto dict_data = read_text_file(dict);
    require(dict_data.ok(), "read sequence dictionary");
    require(dict_data.value().find("@SQ\tSN:chr2\tLN:200\tAS:asm1\tSP:human\tUR:file://ref.fa") != std::string::npos,
            "dictionary first contig with metadata");
    require(dict_data.value().find("chr2") < dict_data.value().find("chr1"), "dictionary preserves fai order");

    auto sam = root / "rg_input.sam";
    write_file(sam,
               "@HD\tVN:1.6\tSO:unknown\n"
               "@RG\tID:old\tSM:oldSample\n"
               "@SQ\tSN:chr1\tLN:1000\n"
               "read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:old\n"
               "read2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
               "read3\t1024\tchr1\t5\t20\t4M\t*\t0\t0\tTGCA\t!!!!\n");
    auto rg_output = root / "rg_output.sam";
    ReadGroupOptions read_group;
    read_group.id = "rg1";
    read_group.sample = "sample1";
    read_group.library = "lib1";
    read_group.platform = "ILLUMINA";
    read_group.platform_unit = "unit1";
    auto rg_result = replace_read_groups_sam_text(sam, rg_output, read_group, true);
    require(rg_result.ok(), "replace read groups");
    auto rg_data = read_text_file(rg_output);
    require(rg_data.ok(), "read read-group output");
    require(rg_data.value().find("@RG\tID:rg1\tSM:sample1\tLB:lib1\tPL:ILLUMINA\tPU:unit1") != std::string::npos,
            "new read-group header");
    require(rg_data.value().find("ID:old") == std::string::npos, "old read-group header removed");
    require(rg_data.value().find("RG:Z:old") == std::string::npos, "old read-group tag removed");
    require(rg_data.value().find("RG:Z:rg1") != std::string::npos, "new read-group tag added");

    auto summary = collect_alignment_summary_text(rg_output);
    require(summary.ok(), "collect alignment summary");
    require(summary.value().total_records == 3, "alignment summary total records");
    require(summary.value().mapped_records == 2, "alignment summary mapped records");
    require(summary.value().unmapped_records == 1, "alignment summary unmapped records");
    require(summary.value().duplicate_records == 1, "alignment summary duplicate records");
    require(summary.value().total_bases == 8, "alignment summary total bases");
    require(summary.value().aligned_bases == 8, "alignment summary aligned bases");
    auto json_summary = format_alignment_summary(summary.value(), true);
    require(json_summary.find("\"mean_mapq\": 40") != std::string::npos, "alignment summary mean mapq");
}

void test_phase3_clean_reorder_fixmate() {
    auto root = temp_root();
    auto fai = root / "ref.fai";
    write_file(fai,
               "chr2\t200\t0\t80\t81\n"
               "chr1\t10\t210\t80\t81\n");
    auto dict = read_fai(fai);
    require(dict.ok(), "read clean/reorder fai");

    auto dirty = root / "dirty.sam";
    write_file(dirty,
               "@HD\tVN:1.6\tSO:unknown\n"
               "@SQ\tSN:chr1\tLN:10\n"
               "@SQ\tSN:chr2\tLN:200\n"
               "good\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n"
               "overflow\t0\tchr1\t9\t60\t5M\t*\t0\t0\tACGTG\t!!!!!\n"
               "missing\t0\tchrMissing\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n");
    auto cleaned = root / "cleaned.sam";
    auto clean_result = clean_sam_text(dirty, cleaned, &dict.value(), true);
    require(clean_result.ok(), "clean sam text");
    auto clean_data = read_text_file(cleaned);
    require(clean_data.ok(), "read cleaned sam");
    require(clean_data.value().find("good\t0\tchr1\t1\t60\t4M") != std::string::npos, "clean preserves good record");
    require(clean_data.value().find("overflow\t4\t*\t0\t0\t*") != std::string::npos, "clean unmaps overflow record");
    require(clean_data.value().find("BC:Z:reference_overflow") != std::string::npos, "clean overflow reason");
    require(clean_data.value().find("BC:Z:reference_conflict") != std::string::npos, "clean missing contig reason");

    auto reordered = root / "reordered.sam";
    auto reorder_result = reorder_sam_text(dirty, reordered, fai, true, true);
    require(reorder_result.ok(), "reorder sam text");
    auto reorder_data = read_text_file(reordered);
    require(reorder_data.ok(), "read reordered sam");
    require(reorder_data.value().find("@SQ\tSN:chr2") < reorder_data.value().find("@SQ\tSN:chr1"), "reorder uses fai order");
    require(!reorder_sam_text(dirty, root / "reorder_fail.sam", fai, false, true).ok(), "reorder rejects missing contig by default");

    auto mate_input = root / "mates.sam";
    write_file(mate_input,
               "@HD\tVN:1.6\tSO:queryname\n"
               "@SQ\tSN:chr1\tLN:100\n"
               "pair1\t65\tchr1\t10\t50\t5M\t*\t0\t0\tAAAAA\t!!!!!\n"
               "pair1\t129\tchr1\t40\t60\t5M\t*\t0\t0\tCCCCC\t!!!!!\n"
               "single\t0\tchr1\t70\t20\t5M\t*\t0\t0\tGGGGG\t!!!!!\n");
    auto fixed = root / "fixed.sam";
    auto fix_result = fix_mate_sam_text(mate_input, fixed, true);
    require(fix_result.ok(), "fix mate sam text");
    auto fixed_data = read_text_file(fixed);
    require(fixed_data.ok(), "read fixed mate sam");
    require(fixed_data.value().find("pair1\t65\tchr1\t10\t50\t5M\t=\t40\t35") != std::string::npos, "fix mate first fields");
    require(fixed_data.value().find("pair1\t129\tchr1\t40\t60\t5M\t=\t10\t-35") != std::string::npos, "fix mate second fields");
    require(fixed_data.value().find("MC:Z:5M\tMQ:i:60") != std::string::npos, "fix mate first tags");
    require(fixed_data.value().find("MC:Z:5M\tMQ:i:50") != std::string::npos, "fix mate second tags");
    require(fixed_data.value().find("single\t0\tchr1\t70") != std::string::npos, "fix mate preserves singleton");
}

void test_phase3_mark_duplicates() {
    auto root = temp_root();
    auto input = root / "dups.sam";
    write_file(input,
               "@HD\tVN:1.6\tSO:queryname\n"
               "@RG\tID:rg1\tLB:lib1\n"
               "@SQ\tSN:chr1\tLN:1000\n"
               "pairKeep\t65\tchr1\t10\t60\t5M\t=\t40\t35\tAAAAA\tIIIII\tRG:Z:rg1\n"
               "pairKeep\t129\tchr1\t40\t60\t5M\t=\t10\t-35\tCCCCC\tIIIII\tRG:Z:rg1\n"
               "pairDup\t65\tchr1\t10\t20\t5M\t=\t40\t35\tAAAAA\t!!!!!\tRG:Z:rg1\n"
               "pairDup\t129\tchr1\t40\t20\t5M\t=\t10\t-35\tCCCCC\t!!!!!\tRG:Z:rg1\n"
               "singleKeep\t0\tchr1\t70\t50\t5M\t*\t0\t0\tGGGGG\tIIIII\tRG:Z:rg1\n"
               "singleDup\t0\tchr1\t70\t20\t5M\t*\t0\t0\tTTTTT\t!!!!!\tRG:Z:rg1\n"
               "secondary\t256\tchr1\t70\t20\t5M\t*\t0\t0\tNNNNN\t!!!!!\tRG:Z:rg1\n"
               "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n");

    auto marked = root / "marked.sam";
    auto metrics_file = root / "mark.metrics.json";
    auto result = mark_duplicates_sam_text(input, marked, &metrics_file, false, true, true);
    require(result.ok(), "mark duplicates sam text");
    require(result.value().records == 8, "mark duplicate records counted");
    require(result.value().eligible_records == 6, "mark duplicate eligible records counted");
    require(result.value().duplicate_sets == 2, "mark duplicate set count");
    require(result.value().duplicates_marked == 3, "mark duplicate marked count");
    auto marked_data = read_text_file(marked);
    require(marked_data.ok(), "read marked duplicates");
    require(marked_data.value().find("pairKeep\t65\tchr1\t10") != std::string::npos, "mark duplicate keeps best pair first end");
    require(marked_data.value().find("pairDup\t1089\tchr1\t10") != std::string::npos, "mark duplicate flags pair first end");
    require(marked_data.value().find("pairDup\t1153\tchr1\t40") != std::string::npos, "mark duplicate flags pair second end");
    require(marked_data.value().find("singleDup\t1024\tchr1\t70") != std::string::npos, "mark duplicate flags singleton");
    require(marked_data.value().find("DT:Z:LB") != std::string::npos, "mark duplicate type tag");
    auto metrics = read_text_file(metrics_file);
    require(metrics.ok(), "read mark duplicate metrics");
    require(metrics.value().find("\"duplicates_marked\": 3") != std::string::npos, "mark duplicate metrics count");

    auto removed = root / "removed.sam";
    auto removed_metrics = root / "removed.metrics.json";
    auto remove_result = mark_duplicates_sam_text(input, removed, &removed_metrics, true, false, true);
    require(remove_result.ok(), "remove duplicates sam text");
    require(remove_result.value().removed_duplicates == 3, "remove duplicate count");
    auto removed_data = read_text_file(removed);
    require(removed_data.ok(), "read removed duplicate output");
    require(removed_data.value().find("pairDup") == std::string::npos, "remove duplicate pair records");
    require(removed_data.value().find("singleDup") == std::string::npos, "remove duplicate singleton");
    require(removed_data.value().find("secondary\t256") != std::string::npos, "remove duplicate preserves secondary");
}

void test_atomic_output_force() {
    auto root = temp_root();
    auto path = root / "out.txt";
    auto first = write_text_file_atomic(path, "one", false);
    require(first.ok(), "first atomic write");
    auto blocked = write_text_file_atomic(path, "two", false);
    require(!blocked.ok(), "blocked overwrite without force");
    auto replaced = write_text_file_atomic(path, "two", true);
    require(replaced.ok(), "overwrite with force");
    auto data = read_text_file(path);
    require(data.ok() && data.value() == "two", "atomic replacement content");
}

} // namespace

int main() {
    test_cigar();
    test_hts_status_and_binary_guards();
    test_chain_index_and_interval_mapping();
    test_reverse_chain_mapping();
    test_fai_validation();
    test_sam_validate_sort_lift();
    test_vcf_lift();
    test_sample_index_and_sam_lift();
    test_phase3_sam_utilities();
    test_phase3_clean_reorder_fixmate();
    test_phase3_mark_duplicates();
    test_atomic_output_force();
    std::cout << "all burnham unit tests passed\n";
    return 0;
}
