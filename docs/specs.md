# Burnham Initial Specs

## Mapping Model

Burnham index files store normalized blocks:

```text
source_contig source_start source_end dest_contig dest_start dest_end strand chain_id
```

Coordinates are 0-based half-open. `strand` is the relationship from source to destination. Position mapping is exact inside one block. Interval mapping requires both interval ends to map to the same destination contig and strand.

## Reason Taxonomy

- `lifted`: record or interval mapped without structural changes.
- `dry_run_lifted`: dry-run equivalent of `lifted`.
- `missing_mapping`: source contig or position is absent from the chain index.
- `split_mapping`: interval spans multiple blocks or maps discontinuously.
- `invalid_interval`: malformed interval coordinates.
- `invalid_cigar`: malformed CIGAR or unsupported validation state.
- `reference_conflict`: record references a contig absent from the provided dictionary.
- `unsupported_format`: command received a binary format before htslib integration.

## Chain Direction

`index-chain` requires `--source-side query` or `--source-side target`. The raw UCSC chain target/query fields are converted into normalized blocks at index-build time. FAI dictionaries are optional in this slice but, when provided, contig names and lengths are checked against the chosen source/destination sides.

## SAM Text Policy

The current SAM tools support text SAM only. `validate-sam` checks headers, reference dictionary consistency, POS, CIGAR syntax, query-length agreement, and reference-bound spans when a dictionary is available. `lift-sam` maps POS and reference names when an alignment start maps; CIGAR transformation, mate repair, and binary formats require the htslib-backed phase. `explain-read` reports the first matching QNAME's source position, mapped destination, chain id, strand, and failure reason.

## Sample-Aware Index Policy

`index-sample` builds a Burnham sample index v1 from text VCF records. The current implementation includes biallelic non-symbolic variants only, optionally filtered by `--sample` and `--haplotype`. Insertions and deletions create deterministic coordinate shifts for downstream SAM text records. Reference positions removed by a selected deletion are rejected with `deleted_by_variant`. Full levioSAM-style haplotype-aware CIGAR transformation is reserved for the htslib-backed phase.

## Phase 3 Text Utility Policy

`dict` creates a SAM sequence dictionary from an FAI file while preserving FAI contig order and optional assembly/species/URI metadata. `read-groups` replaces existing `@RG` headers and `RG:Z` tags with one explicit read group. `alignment-summary` reports deterministic text or JSON counters for records, mapped/unmapped status, secondary/supplementary/duplicate flags, bases, percent mapped, and mean MAPQ. `clean` normalizes existing unmapped records and converts invalid mapped records to unmapped records with `BC:Z:<reason>` tags. `reorder` rewrites `@SQ` headers in FAI order and rejects missing contigs unless `--allow-missing-contigs` is set. `fix-mate` repairs RNEXT/PNEXT/TLEN and MC/MQ tags for simple two-record queryname pairs. `mark-dup` groups eligible primary mapped text records by library and unclipped 5-prime geometry, keeps the highest base-quality score representative with MAPQ and input order as tie-breakers, marks non-representatives with `0x400`, and can emit `DT:Z:LB`, remove duplicate records, and write JSON metrics. These are text SAM foundations for Picard-style behavior; binary parity requires htslib integration.

## VCF Text Policy

The current `lift-vcf-chain` maps POS and CHROM for text VCF records. Allele normalization, left-alignment, swapped ref/alt recovery, genotype remapping, and INFO-field policy are specified as required next work before Picard-compatible VCF liftover is declared complete.
