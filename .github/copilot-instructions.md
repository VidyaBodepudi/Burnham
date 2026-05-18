# Burnham Workspace Instructions

- Use C++20 and CMake for all production code.
- Keep the dependency surface small until htslib/vcpkg integration is added.
- Prefer focused, tested increments over broad unverified rewrites.
- Commands should keep compatibility with the plan in `docs/roadmap.md`.
- Run `cmake --build` and `ctest` after code changes when the local toolchain is available.

## Setup Progress

- [x] Clarify project requirements: Burnham C++20 genomic liftover/post-processing CLI from the supplied plan.
- [x] Scaffold the project structure.
- [x] Implement the first verifiable liftover/platform slice.
- [x] Add tests and documentation.
- [x] Compile and run tests locally with direct LLVM Clang fallback because CMake was unavailable in PATH.
