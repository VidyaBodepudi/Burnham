# Burnham Roadmap

This roadmap condenses the supplied plan into implementation stages.

## Phase 0: Specification and Product Contract

- Normalize mappings into explicit source and destination references.
- Reject ambiguous chain direction unless the user provides `--source-side`.
- Keep a stable reason taxonomy for liftover and validation outcomes.
- Specify mate policy, CIGAR validation, VCF allele/genotype policy, and Picard/levioSAM compatibility targets before parity work.

## Phase 1: Core Platform

- CMake 3.25+, C++20, CMake presets, CI-ready tests.
- Minimal dependency surface first; htslib/vcpkg integration boundary is now scaffolded with `vcpkg.json`, `BURNHAM_ENABLE_HTSLIB`, `vcpkg-debug`, and `doctor` status reporting.
- Atomic output writer, FAI dictionary parser, Result/Error model, metrics counters, and deterministic behavior.
- Bounded streaming/concurrency and full binary format wrappers are next platform milestones.

## Phase 2: Liftover MVP

- Implemented in this slice: `index-chain`, `index-sample`, `inspect-index`, `inspect-sample-index`, `lift-bed`, simple text `lift-sam`, `explain-read`, and simple text `lift-vcf-chain`.
- Implemented for v1 boundary: htslib-backed `.bam`/`.cram` validation probe through `validate-sam` when built with `BURNHAM_ENABLE_HTSLIB=ON`; binary transformation commands still reject `.bam`/`.cram` with explicit text-boundary errors.
- Next: full htslib-backed SAM/BAM/CRAM command ports, KSW2 realignment, mate-pair policy, richer sample haplotype CIGAR transforms, and Picard-style VCF allele normalization.

## Phase 3: Picard-Style Tools

- Implemented in this slice: text SAM `sort`, `validate-sam`, `dict`, `read-groups`, `alignment-summary`, `clean`, `fix-mate`, `reorder`, and `mark-dup` foundations.
- MarkDuplicates algorithm spec is pinned in `docs/mark-duplicates-spec.md`; the current implementation covers library-aware text duplicate sets, quality/MAPQ/input-order representative selection, duplicate flags, optional `DT:Z:LB` tags, removal mode, and JSON metrics.
- Next: binary format support, optical duplicate detection, UMI-aware grouping, platform presets, and Picard parity fixtures.

## Phase 4: Performance Work

- Implemented in this slice: Phase 4 planning document, deterministic text-foundation benchmark harness in `scripts/bench-phase4.ps1`, and median-time regression comparison in `scripts/compare-phase4.ps1`.
- Current text-foundation Phase 4 status: complete. The harness now records repeated samples, median/mean/min/max timings, throughput, input/output sizes, and benchmark directory size.
- Next: baseline levioSAM2/Picard profiles, htslib threading, rank/select benchmarks, temp I/O, allocator profiling, and measured claims after production IO/parity work lands.
- See `docs/phase4-plan.md`, `docs/phase4-completion.md`, and `docs/phase-0-3-audit.md`; production-scale benchmarking still depends on htslib and golden datasets.

## Phase 5: Verification, CI, Packaging

- Implemented in this slice: unit tests for chain parsing/indexing, interval mapping, CIGAR validation, SAM validation/sorting/liftover, VCF liftover, Phase 3 SAM utilities, mark-dup behavior, atomic writes, Phase 5 local verification, GitHub Actions CI scaffolding, CMake install rules, and Windows direct-Clang package creation.
- Current text-foundation Phase 5 status: complete. See `docs/phase5-plan.md` and `docs/phase5-completion.md`.
- Next: htslib-backed wrapper tests, fuzzing, golden parity tests, Docker/conda packaging, migration docs, and signed/provenance release artifacts after production IO/parity work lands.
