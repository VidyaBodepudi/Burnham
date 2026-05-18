# Phase 4 Plan: Profile-Driven Performance

Phase 4 turns Burnham from a verified implementation into a measured implementation. The rule for this phase is simple: every performance claim needs a reproducible command, dataset definition, metric, and baseline.

## Scope

- Measure current text-format commands with repeatable local fixtures.
- Establish the result schema and benchmark workflow before larger dependencies are added.
- Add htslib-backed benchmarks once BAM/CRAM/VCF-BCF support is implemented.
- Compare production paths against pinned Picard, levioSAM, and levioSAM2 fixtures before claiming parity or speedups.

## Completed Text-Foundation Slice

The completed Phase 4 text slice is a benchmark and regression workflow built from two scripts:

- `scripts/bench-phase4.ps1` generates deterministic SAM, BED, and chain fixtures, runs representative Burnham commands, records repeated elapsed-time samples, calculates summary statistics, and writes schema-versioned JSON.
- `scripts/compare-phase4.ps1` compares two benchmark JSON files by command name and flags regressions above a configurable percentage threshold.

This workflow intentionally measures only the current text implementation. It is useful for regression tracking and workflow design, but it is not a substitute for future htslib/BAM/CRAM profiling.

## Benchmark Groups

- Liftover: `index-chain`, `lift-bed`, and text `lift-sam`.
- SAM utilities: `validate-sam`, `sort`, `alignment-summary`, and `mark-dup`.
- Future htslib group: BAM/CRAM parse/decompression, output compression, threaded writing, and binary validation.
- Future parity group: Picard, levioSAM, and levioSAM2 golden comparisons.

## Metrics Schema

Phase 4 benchmark JSON uses schema version 2:

- `schema_version`: integer schema version.
- `phase`: benchmark family, currently `phase4-text-foundation`.
- `generated_at_utc`: ISO-8601 timestamp.
- `record_count`: requested SAM record count for generated fixtures.
- `repeat_count`: measured repetitions per command.
- `warmup_count`: warmup repetitions per command.
- `benchmark_dir`: directory containing generated fixtures and command outputs.
- `benchmark_dir_bytes`: total bytes in the benchmark directory after the run.
- `note`: caveat text for the benchmark slice.
- `commands`: array of command measurements with `name`, `command`, `samples_milliseconds`, `min_milliseconds`, `median_milliseconds`, `mean_milliseconds`, `max_milliseconds`, `records_per_second`, `input_bytes`, `output_bytes`, and `exit_code`.

## Regression Workflow

Create a baseline:

```powershell
./scripts/bench-phase4.ps1 -RecordCount 2000 -RepeatCount 5 -OutputPath build/bench/phase4-baseline.json
```

Create a current run after code changes:

```powershell
./scripts/bench-phase4.ps1 -RecordCount 2000 -RepeatCount 5 -OutputPath build/bench/phase4-current.json -SkipBuild
```

Compare them:

```powershell
./scripts/compare-phase4.ps1 -BaselinePath build/bench/phase4-baseline.json -CurrentPath build/bench/phase4-current.json -RegressionThresholdPercent 15 -FailOnRegression
```

## Development Backlog

1. Add htslib/vcpkg and binary record wrappers.
2. Extend the benchmark harness with BAM/CRAM inputs and compression/threading variants.
3. Add pinned Picard/levioSAM/levioSAM2 baseline runners when those tools and fixtures are available locally.
4. Add rank/select and CIGAR-transform microbenchmarks once the production engines exist.
5. Add temp-space, peak-memory, and per-stage throughput metrics.
6. Gate performance claims on reproducible benchmark JSON committed or attached with release artifacts.

## Exit Criteria

- Local text benchmark harness runs after a clean build.
- Text benchmark JSON includes repeat samples, median/mean/min/max timings, throughput, and input/output sizes.
- Local regression comparator can fail builds when median time regresses past an explicit threshold.
- Binary IO benchmark group exists after htslib integration.
- Golden baseline commands are documented with pinned tool versions and fixture sources.
- Reported speedups include wall time, CPU time where available, memory/temp-space notes, dataset definitions, and command lines.