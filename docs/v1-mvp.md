# V1 MVP: Production IO Boundary

The v1 MVP builds on the completed text-format foundation and adds the first production IO boundary needed for htslib-backed development.

## Completed in This Slice

- `vcpkg.json` declares `htslib` as the production IO dependency.
- `CMakePresets.json` includes `vcpkg-debug` for htslib-enabled local configuration.
- `CMakeLists.txt` exposes `BURNHAM_ENABLE_HTSLIB` and links a known htslib CMake target when enabled.
- `burnham doctor` reports htslib build status.
- `validate-sam` can route `.bam` and `.cram` inputs through an htslib read probe when Burnham is built with htslib enabled.
- Text-only SAM transformation commands now reject `.bam` and `.cram` inputs with an explicit boundary message instead of treating binary data as malformed text.
- `scripts/verify-v1-mvp.ps1` runs the local MVP gate.

## Local Test Command

```powershell
powershell -ExecutionPolicy Bypass -File scripts\verify-v1-mvp.ps1
```

The summary is written to `build/v1-mvp/verification.json`.

## htslib-Enabled Build Path

When CMake and vcpkg are available, configure with:

```powershell
cmake --preset vcpkg-debug
cmake --build --preset vcpkg-debug
ctest --preset vcpkg-debug
```

This path is expected to become the main production build after htslib-backed command ports land.

## Remaining v1 Work

- Port `lift-sam`, `sort`, `mark-dup`, `fix-mate`, `alignment-summary`, and validation internals from text records to shared htslib-backed alignment record adapters.
- Add real BAM/CRAM fixtures and golden parity comparisons against Picard, levioSAM, and levioSAM2.
- Add CRAM reference handling policy.
- Extend Phase 4 benchmarks to binary parse/decompression, compression/threading, and htslib writer paths.