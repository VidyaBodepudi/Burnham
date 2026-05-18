# Phase 0-3 Audit

This audit captures the state after completing the Phase 3 text foundation and before starting Phase 4 performance work.

## Status

- Phase 0 specification: complete for the current implementation slice. Mapping direction, coordinate normalization, reason taxonomy, SAM text policy, sample-index behavior, and MarkDuplicates foundation behavior are documented.
- Phase 1 platform: complete for local text development. The project has a C++20/CMake scaffold, a small dependency surface, Result/Error handling, atomic output, FAI parsing, deterministic text IO, and a verified direct LLVM Clang build fallback.
- Phase 2 liftover MVP: complete for text inputs. Implemented commands include `index-chain`, `index-sample`, `inspect-index`, `inspect-sample-index`, `lift-bed`, text `lift-sam`, `explain-read`, and text `lift-vcf-chain`.
- Phase 3 Picard-style tools: complete for text SAM foundations. Implemented commands include `sort`, `validate-sam`, `dict`, `read-groups`, `alignment-summary`, `clean`, `reorder`, `fix-mate`, and `mark-dup`.

## Remaining Before Production Parity

- htslib/vcpkg integration for BAM/CRAM/VCF-BCF IO and BGZF threading.
- Golden parity fixtures against Picard and levioSAM behavior.
- Full mate-pair liftover policy and haplotype-aware CIGAR/reference transformation.
- KSW2 or equivalent realignment support for indel/variant-aware alignment repair.
- Full VCF allele normalization, left-alignment, swapped REF/ALT recovery, genotype remapping, and INFO policy.
- MarkDuplicates optical duplicate detection, UMI-aware grouping, platform presets, and full compatibility-profile validation.
- CI coverage with sanitizers, fuzzing, packaging checks, and real CMake availability in the local/dev environment.

## Phase 4 Readiness

The text foundation is ready to move into Phase 4-style measurement for local development: commands are deterministic, tested, and runnable with `scripts/test-clang.ps1`. For production-scale Phase 4 claims, integrate htslib and golden datasets first so profiling measures the real binary IO, threading, and parity surfaces instead of only the current text slice.