# MarkDuplicates Specification

This document pins Burnham's duplicate-marking contract before implementation. It is intentionally written as an algorithm spec, not as a port of Picard internals.

## Scope

`mark-dup` marks duplicate SAM/BAM/CRAM records by setting flag `0x400` and optional duplicate-type tags. The first implementation may operate on text SAM, but the behavioral contract must be compatible with future htslib-backed records.

Out of scope for the first implementation:

- Cohort-level duplicate marking across files.
- Flow-based and Ultima-specific scoring beyond an explicit platform preset placeholder.
- Optical duplicate image geometry beyond the parser and distance contract below.
- UMI-aware duplicate families until UMI tags and error rules are pinned.

## Eligibility

Records are ineligible for duplicate marking when any of these are true:

- Unmapped flag `0x4` is set.
- Secondary flag `0x100` is set.
- Supplementary flag `0x800` is set, unless a later compatibility profile explicitly opts in.
- The record has malformed required SAM fields.
- The record has an invalid CIGAR.

Ineligible records pass through unchanged except for validation errors in strict mode.

## Duplicate Keys

Burnham groups eligible records by library and alignment geometry.

Library key:

- Prefer `RG:Z` tag mapped through the `@RG` header's `LB` value.
- If `RG:Z` is present but no matching `@RG`/`LB` exists, use the read group ID.
- If no `RG:Z` exists, use the literal library key `unknown`.

Single-end key:

- Library key.
- Reference name.
- 5-prime unclipped position.
- Read strand.

Paired-end key:

- Library key.
- First end reference name and 5-prime unclipped position.
- First end strand.
- Mate reference name and mate 5-prime unclipped position when available.
- Mate strand when available.

For the text SAM foundation, mate 5-prime position may use `PNEXT` and mate strand flag `0x20`. The htslib implementation must compute the mate geometry from actual mate records when both ends are present.

## 5-Prime Position

The unclipped 5-prime position is derived from POS and CIGAR:

- Forward strand: subtract leading hard/soft clips from POS.
- Reverse strand: add reference span and trailing hard/soft clips, then subtract one to get the 5-prime coordinate.
- Coordinates remain 1-based for SAM-facing reporting and deterministic key construction.

If CIGAR is `*`, the record is ineligible.

## Representative Selection

Within each duplicate set, Burnham keeps one representative unmarked.

Scoring order:

1. Higher sum of base qualities after subtracting ASCII 33, ignoring bases with quality below 15.
2. Higher MAPQ.
3. Earlier input order.

All non-representative eligible records in the set are marked duplicate with flag `0x400`.

## Optical Duplicate Contract

The first implementation may mark all duplicates as library duplicates. Optical duplicate classification is enabled only when a read-name parser is configured.

Default read-name parser:

- Illumina-style names split on `:` with lane/tile/x/y in the last four numeric fields.
- Optical duplicates share library key, duplicate key, lane, and tile.
- Pixel distance is `sqrt((x1-x2)^2 + (y1-y2)^2)`.
- Default optical distance threshold is 100 for Illumina-style short reads.

Optical duplicate records receive duplicate type tag `DT:Z:SQ`; library duplicates receive `DT:Z:LB` when duplicate type tagging is enabled.

## UMI Contract

UMI-aware duplicate marking is disabled by default. When added, UMI mode must specify:

- Accepted UMI tag names.
- Duplex versus single UMI behavior.
- Error tolerance and normalization.
- Whether missing UMI records are grouped together or treated as ineligible.

No UMI behavior is implicit in the first implementation.

## Platform Presets

Initial accepted preset values are `illumina`, `ultima`, `pacbio`, and `ont`.

For the first implementation:

- `illumina`: enables Illumina optical parser defaults when optical classification is requested.
- `ultima`, `pacbio`, `ont`: do not enable optical duplicate classification by default.

Platform-specific scoring changes require separate benchmarks and parity fixtures.

## Output

Default behavior:

- Preserve input order.
- Preserve all non-duplicate flags and tags.
- Set `0x400` on duplicate records.
- Do not remove duplicates unless `--remove-duplicates` is set.
- Emit JSON metrics when `--metrics-file` is provided.

Metrics schema v1 fields:

- `records`
- `eligible_records`
- `duplicate_sets`
- `duplicates_marked`
- `library_duplicates`
- `optical_duplicates`
- `removed_duplicates`

## Compatibility Profiles

`--compat picard` aims for semantic parity with the pinned Picard version once golden fixtures exist. `--compat burnham-strict` may reject malformed records earlier and must document every intentional difference.

## Verification Gates

- Unit tests for CIGAR unclipping, duplicate keys, base-quality scoring, tie-breaks, flag/tag updates, and metrics.
- Golden tests against pinned Picard for representative short-read paired-end datasets.
- Edge cases for unmapped, secondary, supplementary, duplicate-preflagged, missing read groups, invalid CIGARs, reverse-strand reads, equal-score tie-breaks, and duplicate removal.
- Deterministic output across repeated runs and different thread counts.
