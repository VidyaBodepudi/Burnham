# Phase 5 Plan: Verification, CI, Packaging

Phase 5 makes Burnham repeatable outside a single developer session. The current slice focuses on the verified text-format foundation while leaving production htslib/golden-fixture expansion explicit.

## Scope

- Keep local verification one command away.
- Add CI scaffolding for CMake builds, tests, sanitizer checks, and the Windows direct-Clang fallback.
- Package the current Windows text-foundation executable with docs and checksums.
- Preserve a clear backlog for fuzzing, golden parity fixtures, binary-format validation, and release distribution.

## Current Development Slice

- `scripts/verify-phase5.ps1` runs unit tests, Phase 4 benchmark smoke checks, benchmark comparison smoke checks, and packaging smoke checks.
- `scripts/package-release.ps1` creates a Windows x64 package from the direct-Clang build, including `burnham.exe`, README, docs, checksums, a zip archive, and a manifest.
- `.github/workflows/ci.yml` defines CMake builds on Linux, macOS, and Windows; ASan/UBSan checks on Ubuntu Clang; and the Windows direct-Clang Phase 5 verification path.

## Local Verification

```powershell
powershell -ExecutionPolicy Bypass -File scripts\verify-phase5.ps1
```

The summary is written to `build/phase5/verification.json`.

## Packaging

```powershell
powershell -ExecutionPolicy Bypass -File scripts\package-release.ps1
```

The default package directory is `build/package/`.

## Backlog Before Production Release

- Add htslib-backed BAM/CRAM/VCF wrapper tests.
- Add fuzz targets for chain, SAM headers, VCF headers, CIGAR strings, and aux tags.
- Add pinned Picard, levioSAM, and levioSAM2 golden fixtures.
- Add clang-tidy and clang-format enforcement once formatting policy is pinned.
- Add Docker and conda/bioconda packaging when binary IO dependencies are available.
- Add release signing/provenance after package contents stabilize.

## Exit Criteria

- Local Phase 5 verification passes.
- Package script produces a zip, manifest, docs, and checksum file.
- CI workflow captures cross-platform CMake, sanitizer, and direct-Clang fallback paths.
- Remaining production verification gaps are documented rather than hidden.