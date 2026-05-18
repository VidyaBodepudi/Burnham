# Burnham

Burnham is a C++20 genomic alignment post-processing and liftover command-line tool. This repository currently contains verified implementation slices from the supplied product plan: a small dependency-light platform, chain indexing, sample-aware VCF indexing, index inspection, BED liftover, text SAM validation/sorting/liftover, simple VCF coordinate liftover, read-group replacement, sequence dictionary creation, alignment summaries, duplicate marking, explain-read, JSON-like metrics, and tests.

The htslib-backed SAM/BAM/CRAM and full Picard/levioSAM parity work is staged in `docs/roadmap.md` and `docs/specs.md`. The current code is intentionally strict about its capabilities rather than pretending to support binary genomics formats before htslib is integrated.

Duplicate marking is implemented for text SAM according to the foundation contract in `docs/mark-duplicates-spec.md`. Optical duplicate detection, UMI-aware grouping, and binary BAM/CRAM parity remain staged for the htslib-backed phase.

## Build

```powershell
cmake --preset default
cmake --build --preset default
ctest --preset default
```

If Ninja is not installed, configure with another generator:

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

On this Windows workspace, CMake was not available in PATH during initial verification, but LLVM Clang was installed. The direct-build fallback is:

```powershell
./scripts/test-clang.ps1
```

## Phase 4 Benchmarks

The current Phase 4 harness measures the deterministic text-format foundation and writes schema-versioned JSON. It is useful for local regression tracking while htslib/BAM/CRAM performance work is still pending.

```powershell
./scripts/bench-phase4.ps1 -RecordCount 2000 -RepeatCount 5
```

The default output is `build/bench/phase4/phase4-benchmark.json`.

Compare a baseline and current run with:

```powershell
./scripts/compare-phase4.ps1 -BaselinePath build/bench/phase4-baseline.json -CurrentPath build/bench/phase4-current.json -RegressionThresholdPercent 15 -FailOnRegression
```

## Phase 5 Verification and Packaging

Run the current text-foundation verification gate with:

```powershell
./scripts/verify-phase5.ps1
```

Create a Windows x64 direct-Clang package with:

```powershell
./scripts/package-release.ps1
```

The package script writes a zip archive, manifest, docs, and SHA-256 checksum under `build/package/`.

## V1 MVP Gate

The current v1 MVP gate verifies the completed text foundation plus the htslib/vcpkg production IO boundary:

```powershell
./scripts/verify-v1-mvp.ps1
```

For an htslib-enabled CMake build, set `VCPKG_ROOT` and use:

```powershell
cmake --preset vcpkg-debug
cmake --build --preset vcpkg-debug
ctest --preset vcpkg-debug
```

## Commands

```text
burnham doctor [--chain-index index.bci] [--sample-index index.bsi]
burnham index-chain --chain input.chain --output output.bci --source-side query|target [--source-fai source.fai] [--dest-fai dest.fai]
burnham index-sample --vcf input.vcf --output output.bsi [--sample SAMPLE] [--haplotype all|0|1]
burnham inspect-index --index output.bci [--json]
burnham inspect-sample-index --index output.bsi [--json]
burnham lift-bed --index output.bci --input input.bed --output lifted.bed [--unmapped-output rejected.bed] [--metrics-file metrics.json] [--dry-run] [--force]
burnham validate-sam --input input.sam [--reference-fai ref.fai] [--summary]
burnham sort --input input.sam --output sorted.sam --order coordinate|queryname [--force]
burnham dict --reference-fai ref.fai --output ref.dict [--assembly ASSEMBLY] [--species SPECIES] [--uri URI] [--force]
burnham read-groups --input input.sam --output output.sam --id ID --sample SAMPLE [--library LIB] [--platform PLATFORM] [--platform-unit UNIT] [--force]
burnham alignment-summary --input input.sam [--output metrics.txt] [--json] [--force]
burnham clean --input input.sam --output clean.sam [--reference-fai ref.fai] [--force]
burnham reorder --input input.sam --output reordered.sam --reference-fai ref.fai [--allow-missing-contigs] [--force]
burnham fix-mate --input input.sam --output fixed.sam [--force]
burnham mark-dup --input input.sam --output marked.sam [--metrics-file metrics.json] [--remove-duplicates] [--duplicate-type-tags] [--force]
burnham lift-sam --index output.bci --input input.sam --output lifted.sam [--unmapped-output rejected.sam] [--reason-tags] [--metrics-file metrics.json] [--dry-run] [--force]
burnham lift-sam --sample-index output.bsi --input input.sam --output lifted.sam [--unmapped-output rejected.sam] [--reason-tags] [--metrics-file metrics.json] [--dry-run] [--force]
burnham explain-read --index output.bci --input input.sam --qname read-name
burnham lift-vcf-chain --index output.bci --input input.vcf --output lifted.vcf [--rejected-output rejected.vcf] [--metrics-file metrics.json] [--dry-run] [--force]
```

## Current Format Support

- Chain input: UCSC chain text files.
- Chain index: Burnham text index v1 (`.bci`).
- Sample index: Burnham text sample index v1 (`.bsi`) built from biallelic text VCF variants and optional sample/haplotype genotype selection.
- BED input/output: BED3+ text, 0-based half-open intervals.
- SAM input/output: SAM text transformation commands in this implementation slice, with chain-index and sample-index coordinate modes plus text sort, validation, cleaning, reordering, mate fixing, duplicate marking, read-group replacement, dictionary, and alignment-summary utilities. Optional htslib builds add `.bam`/`.cram` validation probing through `validate-sam`; full binary transformation command ports are the next v1 milestone.
- VCF input/output: VCF text coordinate liftover only; allele normalization is specified but not yet implemented.

## Safety Defaults

- `index-chain` requires an explicit `--source-side query|target` to avoid guessing chain direction.
- Output files are written through temporary files and atomically renamed on success.
- Existing outputs require `--force`.
- Liftover commands emit per-reason metrics and can run with `--dry-run`.
