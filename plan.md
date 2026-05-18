# Plan: Burnham — Revised Scope, Stability, and User Delight

## TL;DR
Build Burnham as a single C++20 genomic alignment post-processing and liftover tool, but stage it deliberately. The core product should first become a safer, faster, more observable replacement for levioSAM/levioSAM2: chain-based SAM/BAM liftover, sample/VCF-aware SAM/BAM liftover, chain-based VCF liftover, BED liftover, indexing, validation gates, and excellent metrics. Picard-style tools should follow in tiers, with SortSam and ValidateSamFile first, MarkDuplicates only after a written algorithm spec and parity tests.

The revised plan removes avoidable ambiguity, avoids a giant first release, adds correctness invariants, and adds user-delight features such as doctor, dry-run, explain, inspect-index, reason tags, atomic outputs, and clear compatibility profiles.

## Review Findings

1. Chain direction is a correctness risk. The prior plan hard-coded levioSAM2's interpretation of chain target/query fields. Burnham should not rely on memory or convention here. It needs an internal normalized model: source reference, destination reference, raw chain header fields, and explicit source/destination dictionaries. Index build must validate both sides and reject ambiguous direction unless the user explicitly opts into inversion.

2. The first plan was overscoped. Shipping liftover plus 10+ Picard tools as one initial milestone would delay value and create a huge correctness burden. Picard tools should be tiered after liftover, and each complex tool needs an algorithm spec before implementation.

3. VCF liftover was underspecified. Left-alignment, genotype remapping, symbolic alleles, INFO field preservation, swapped ref/alt recovery, and rejected variant reporting need explicit rules before code is written.

4. SAM/BAM liftover needs hard correctness gates. Mate-pair handling, CIGAR validation, reference-bound checks, unmapped/deferred behavior, and per-read failure reasons must be first-class, not afterthoughts.

5. The SDSL dependency is a maintenance risk. Burnham only needs a small rank/select subset. The revised plan prefers a minimal vendored rank/select implementation first, with SDSL retained only as a benchmark/reference option.

6. The concurrency plan needed order and backpressure contracts. A work-stealing pool is useful for CPU stages, but streaming genomics tools also need bounded queues, sequence numbers for order preservation, queue closure semantics, and clear per-subcommand topologies.

7. Observability should be fast by default. JSON and stderr summaries should ship first with lock-free/low-contention counters. Prometheus/OpenTelemetry should be optional later, because tracing in hot per-record paths can become overhead and API churn.

8. Performance claims need measured baselines. The plan should require profiling levioSAM2/Picard on realistic WGS/WES/long-read data before committing to 2x or 3x promises. Optimizations should follow profiles, not taste.

9. Some commands are redundant for MVP. collate, reconcile, duplicate sort order, conversion helpers, and advanced workflow utilities are useful, but they should not block the core liftover release.

10. User delight was underrepresented. Burnham should make dangerous workflows easy to inspect: dry-run, doctor, inspect-index, explain-read, reason tags, compatibility modes, presets, atomic writes, and actionable summaries.

## Revised Architecture

Burnham remains one binary, but not one overloaded command.

Liftover group:
- index-chain — build/load chain index with explicit source/destination validation
- index-sample — build/load sample-aware VCF bitvector index
- lift-sam — SAM/BAM/CRAM liftover using either chain index or sample index
- lift-vcf-chain — Picard-style chain VCF liftover
- lift-bed — BED interval liftover
- inspect-index — report contig coverage, chain direction, index metadata
- explain-read — trace why a read lifted, changed CIGAR, deferred, or became unmapped

Alignment group, staged after liftover:
- sort — coordinate/queryname sorting first; duplicate sort can be advanced later
- validate-sam — SAM/BAM/CRAM validation with summary and verbose modes
- alignment-summary — Picard-style CollectAlignmentSummaryMetrics equivalent; renamed from metrics to avoid confusion with observability metrics
- mark-dup — Picard-compatible duplicate marking after algorithm spec
- merge, reorder, read-groups, dict, clean, fix-mate, sam2fastq, fastq2sam, convert — later utility tiers

Global UX:
- --dry-run — process input and emit metrics without writing final outputs
- --force — required to overwrite existing outputs
- --progress auto|never|always — progress with records/sec and stage bottleneck
- --metrics-file — JSON metrics schema
- --compat picard|leviosam2|burnham-strict — define output behavior and edge-case policy
- --preset ultima|illumina|ont|pacbio — set sensible platform defaults
- --reason-tags — write per-record/per-variant failure or transformation reason tags when requested

## Steps

### Phase 0: Specification and Product Contract
1. Define Burnham's normalized mapping model: source reference, destination reference, raw chain fields, strand handling, and inversion policy. Index build must validate source and destination dictionaries and write this metadata into the index.
2. Define error/reason taxonomy for all liftover outcomes: lifted unchanged, lifted CIGAR changed, lifted reversed, deferred by rule, unmapped by missing contig, unmapped by out-of-bounds, invalid CIGAR, reference conflict, ambiguous chain, rejected VCF reason.
3. Define SAM mate-pair policy for lift-sam: both lifted, one lifted, neither lifted, pre-existing unmapped, secondary/supplementary records, and how flags/RNEXT/PNEXT/TLEN/MAPQ/MC are updated.
4. Define CIGAR validation contract: no zero-length ops, query length matches SEQ, reference span in bounds, unsupported ops handled predictably, extended ops normalized only when compatibility mode requires it.
5. Define VCF liftover policy: Picard/GATK-style left-alignment, reference validation mandatory by default, genotype remapping rules, INFO preservation/drop rules, symbolic allele behavior, rejected VCF annotations, and deterministic multi-allelic handling.
6. Pin compatibility targets: levioSAM2 version/commit for lift-sam parity, levioSAM version/commit for sample-aware liftover, Picard version/commit for lift-vcf-chain/sort/validate/mark-dup parity.
7. Write MarkDuplicates, ValidateSamFile, and alignment-summary algorithm specifications before implementing those tools.

### Phase 1: Core Platform
1. Build system: CMake 3.25+, C++20, CMakePresets, vcpkg manifest with pinned builtin baseline, binary cache strategy, and CI-tested triplets for Linux x86_64, Linux ARM64, macOS ARM64, and Windows/WSL where practical.
2. Dependencies: htslib via vcpkg; CLI11, spdlog/fmt, nlohmann-json, GoogleTest, benchmark. Vendor KSW2 and a minimal rank/select implementation. Avoid SDSL as a required runtime dependency unless benchmarks prove it superior.
3. Core wrappers: owning and non-owning BAM/VCF record wrappers, type-safe aux/tag accessors, reference dictionary model, FASTA/FAI wrapper, temporary-file manager, atomic output writer, and Result/Error model compatible with C++20.
4. Concurrency: bounded Reader -> Processor -> Writer pipelines with backpressure, sequence numbers for order-preserving stages, graceful shutdown, and per-stage queue metrics. Use work-stealing only inside CPU-heavy processor pools.
5. Metrics: low-contention counters/timers, schema-versioned JSON, stderr summary, per-stage throughput, per-contig summaries, and reason breakdowns. Defer OpenTelemetry/Prometheus integration until after the JSON schema stabilizes.
6. External sort: capped memory default, temp-space preflight, deterministic comparators, compressed temp chunks, automatic cleanup, and crash-safe partial-output handling.
7. CLI: consistent option names, aliases for legacy users, shell completions, --dry-run, --compat, --preset, --force, --progress, and clear error messages.
8. Stability features: fuzz parsers for chain/VCF/SAM headers, validate index metadata, atomic writes, tempfile cleanup, deterministic random seeds for tie-breaks, and reproducible output behavior under multi-threading.

### Phase 2: Liftover MVP
1. index-chain and inspect-index: build chain index, normalize source/destination metadata, validate dictionaries, report contig coverage, gaps, strand distribution, missing contigs, and suspected direction mismatch.
2. index-sample: build sample/haplotype-aware VCF bitvector index from levioSAM's VcfMap concept, with explicit sample/haplotype metadata and variant statistics.
3. lift-sam with chain engine: port levioSAM2 ChainMap behavior, add correctness gates, mate policy, CIGAR validation, deferred/unmapped reason tracking, UG tag support, KSW2 realignment modes, and JSON/stderr metrics.
4. lift-sam with sample engine: port levioSAM VcfMap behavior under the same lift-sam UX. The user chooses engine by providing a chain index or sample index; Burnham does not silently mix them.
5. lift-vcf-chain: implement Picard-style VCF liftover using chain mappings, target reference validation, left-alignment, deterministic genotype/INFO policy, rejected VCF output, per-sample rejection metrics, and symbolic allele policy.
6. lift-bed: port BED interval liftover with clear unmapped interval reporting and coverage metrics.
7. explain-read and --dry-run: make liftover debuggable. Users can inspect one read or run a whole file to see how many reads would change, defer, or fail before writing output.
8. Optional but not MVP-blocking: collate and reconcile remain planned advanced liftover utilities, but they do not block Phase 2 release.

### Phase 3: Picard-Style Tools in Tiers
1. Tier 1: sort and validate-sam. These are high-value, deterministic, and establish the external-sort and validation frameworks.
2. Tier 2: dict, read-groups, clean, fix-mate, reorder. These are mostly streaming utilities and improve daily ergonomics.
3. Tier 3: alignment-summary. Implement after validation/reference wrappers are stable so mismatch and chimera metrics are reliable.
4. Tier 4: mark-dup. Implement only after the MarkDuplicates spec is pinned: duplicate scoring strategy, optical duplicate parsing, optical distance threshold, library/read-group behavior, UMI mode, flow/Ultima mode, tie-breaking, output tags, and parity expectations.
5. Tier 5: merge, sam2fastq, fastq2sam, convert, mark-dup-cohort. Useful, but should not delay the core product.

### Phase 4: Performance Work, Profile-Driven
1. Establish baselines against levioSAM2 and Picard on realistic datasets: WGS, WES, long-read, high-CIGAR-complexity, and Ultima-specific BAMs.
2. Profile before optimizing: htslib parse/decompression, CIGAR lifting, rank/select lookup, KSW2 realignment, output compression, external sort temp I/O, queue contention.
3. Optimize confirmed bottlenecks only. Candidate order: htslib thread allocation, batch rank/select lookups, external-sort temp compression, arena allocation for BAM batches, hot hash map replacement, SIMD reverse-complement, io_uring only if Linux I/O is proven limiting.
4. Set performance claims after measurement, with conditions and dataset definitions. Prefer honest targets over marketing numbers.
5. Add resumable/checkpointed liftover only after the batch writer and atomic output model are stable.

### Phase 5: Verification, CI, Packaging
1. Unit tests: chain parsing/direction, rank/select, CIGAR transforms, mate policy, VCF allele transforms, genotype remapping, metrics schema, external sort, tempfile cleanup, BAM/VCF wrappers.
2. Fuzz tests: malformed chain, VCF, SAM headers, CIGAR strings, aux tags, huge contig names, truncated inputs.
3. Golden tests: compare against pinned levioSAM2, levioSAM, and Picard outputs semantically, not only byte-for-byte. Define acceptable tie-break differences explicitly.
4. Genomics edge cases: secondary/supplementary alignments, unmapped mates, multi-mapped reads, reverse-strand chains, tandem-repeat indels, symbolic alleles, multi-allelic variants, mitochondrial contigs, non-human references, long-read CIGARs, polyploid genotypes, circular-reference disclaimer.
5. CI: Linux x86_64 GCC/Clang, Linux ARM64, macOS ARM64 NEON/scalar, Windows/WSL smoke tests, ASan/UBSan/TSan, clang-tidy, clang-format, coverage gates.
6. Packaging: static-ish release binaries where feasible, Docker image with documented mounted volumes, vcpkg binary cache, conda/bioconda packaging later if demand exists.
7. Documentation: migration guide from levioSAM2/Picard, compatibility profiles, index direction guide, VCF liftover caveats, platform presets, metrics schema, and troubleshooting with doctor/explain examples.

## Relevant Code to Reuse or Reimplement

From levioSAM2:
- chain::ChainMap parser/index/lift_pos/lift_cigar/lift_aln should be the starting point for chain liftover, but wrapped in Burnham's normalized mapping model and correctness gates.
- cigar helpers should be ported but extended with validate_cigar and reference-bound checks.
- WriteDeferred, flag/unmap helpers, KSW2 realignment, UG tag handling, BED utilities, collate, and reconcile should be reused where they remain semantically valid.

From levioSAM:
- lift::Lift and lift::LiftMap provide the sample/VCF-aware bitvector model. Burnham should port this as a distinct sample engine rather than blending it with chain liftover.

From Picard:
- LiftoverVcf/LiftoverUtils should be algorithmic references for chain VCF liftover, left-alignment, swapped ref/alt handling, genotype remapping, and rejection reasons.
- SortSam, ValidateSamFile, MarkDuplicates, MergeSamFiles, CollectAlignmentSummaryMetrics, AddOrReplaceReadGroups, ReorderSam, CreateSequenceDictionary should be reimplemented from written specs, not mechanically ported.

## Verification Gates

1. Direction safety: index-chain rejects ambiguous source/destination metadata and passes known coordinate tests against pinned UCSC/Picard/levioSAM2 fixtures.
2. SAM liftover correctness: positions, flags, mate fields, CIGARs, sequences, qualities, tags, and deferred/unmapped outcomes match pinned levioSAM2 behavior or documented Burnham-strict behavior.
3. Sample-aware liftover correctness: sample/haplotype VCF index produces expected coordinate/CIGAR transforms against levioSAM fixtures.
4. VCF chain liftover correctness: positions, alleles, left-alignment, genotype policy, INFO behavior, and rejected variants match pinned Picard behavior where compatible.
5. Output validity: every emitted BAM/VCF from integration tests passes Burnham validation and, where applicable, Picard/htslib validation.
6. Determinism: multi-threaded runs produce semantically identical outputs across repeated runs.
7. Metrics: schema validation, per-reason counters, per-stage timers, per-contig summaries, and dry-run counts match expected fixture values.
8. Performance: report measured wall-clock, CPU, memory, and temp-space use against pinned baselines before claiming speedups.
9. Stability: sanitizer-clean, fuzz-clean for agreed corpus duration, no tempfile leaks, atomic output behavior verified under interrupted runs.

## Decisions Updated

- Keep Burnham as one binary, but use explicit commands and modes instead of hidden algorithm selection.
- Core MVP must include both chain-based and sample/VCF-aware SAM/BAM liftover, because that is the main product differentiator.
- Chain VCF liftover is included in the liftover tier because it closes an important Picard gap and shares chain/reference infrastructure.
- Picard-style alignment tools are staged after liftover, not part of the first all-or-nothing milestone.
- Observability starts with stderr and JSON. Prometheus/OpenTelemetry are later optional exporters.
- Replace required SDSL dependency with minimal rank/select unless benchmarks prove SDSL is worth carrying.
- Use C++20, but do not rely on std::expected; use a Burnham Result type or tl::expected-style compatibility wrapper until C++23 is adopted.
- Rename the Picard metrics command to alignment-summary to avoid confusion with Burnham's observability metrics.
- Keep excluded scope unchanged: no read mapper, no variant caller, no BQSR, no deprecated GATK indel realignment.

## User-Delight Augmentations

1. doctor: validate installation, htslib features, reference/index compatibility, writable temp/output paths, thread settings, and likely chain-direction mistakes.
2. inspect-index: show contig coverage, missing contigs, strand distribution, chain gaps, index source metadata, and command used to build it.
3. explain-read: given a QNAME, show source interval, destination interval, CIGAR transform, mate policy, reason tags, and whether realignment was used.
4. dry-run: run full logic without writing output, producing rejection/defer/change summaries and estimated output sizes.
5. reason tags and reason summaries: make failures debuggable without reading logs.
6. compatibility modes: Picard-compatible, levioSAM2-compatible, and Burnham-strict behavior for users who need reproducibility or stricter safety.
7. platform presets: illumina, ultima, pacbio, ont settings for duplicate marking, tag handling, realignment, and optical duplicate assumptions.
8. atomic outputs: write to temporary files and rename only on success, with --force required for overwrite.
9. progress with stage bottlenecks: records/sec, read/lift/write timings, compression ratio, temp-space use.
10. clear migration docs: Picard and levioSAM2 command equivalents with Burnham commands and behavior differences.

## Redundancies Removed or Deferred

- collate and reconcile move out of the liftover MVP and become advanced utilities.
- duplicate sort order is not required for initial sort; coordinate/queryname come first.
- OpenTelemetry/Prometheus do not block initial observability.
- io_uring, SIMD, and lock-free everything are not assumed wins; they wait for profiling.
- bam-to-cram/cram-to-bam can be a later convert command because most commands should already support output format selection via htslib.
- sam2fastq/fastq2sam remain useful but are not core to Burnham's differentiating liftover story.
