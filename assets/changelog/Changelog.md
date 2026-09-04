# isotools Changelog

## v0.0.42

**BREAKING CHANGE: iso-classify v0.0.13 — generalizes TOGA-specific naming to reference sets; iso-intron v0.0.17 — matching TSV schema**

- Iso-classify v0.0.13 renames the CLI flags of the intron subcommand: `--isoseq`/`-i` is now `--input`/`-i` (paths to input reads BED12, no longer Iso-Seq-specific) and `--toga`/`-t` is now `--reference`/`-r` (path to reference annotation .bed file). Existing command lines using the old long names must switch to `--input` and `--reference`.
- The per-intron schema field `is_toga_supported` is renamed to `is_reference_supported`, and its textual representation changes from `TOGA_SUPPORT`/`NO_TOGA_SUPPORT` to `REFERENCE_SUPPORT`/`NO_REFERENCE_SUPPORT` (doc comments updated accordingly). The field's semantics are unchanged: it still reports whether the intron is supported by the annotation passed via the reference flag.
- Iso-intron v0.0.17 updates its TSV reader to parse the renamed field (`is_reference_supported`, matching `REFERENCE_SUPPORT`/`NO_REFERENCE_SUPPORT`), so both crates agree on the new schema. Parsing the old `TOGA_SUPPORT` tokens will now fail.
- Bumped the workspace to `0.0.42`, `iso-classify` to `0.0.13` and `iso-intron` to `0.0.17` (`Cargo.toml`/`Cargo.lock`).

## v0.0.41

**iso-intron v0.0.16 — `--allow-missing` covers chromosomes missing from the reference index**

- Iso-intron v0.0.16 extends `--allow-missing` to the chromosome-level index lookup. Previously, when a chromosome had intron-bearing reads but was absent from the reference intron index, iso-intron logged an `ERROR` and exited; with the flag set, that case now logs a `WARN` and continues with an empty index, so all reads on that chromosome are reported free of IRs. Without the flag the strict error and exit are preserved.
- The `WARN` for the no-introns edge case now includes the chromosome name (`No introns found for <chr> in input -> ...`), so multi-chromosome runs no longer log identical messages that are impossible to attribute.
- Bumped `iso-intron` to `0.0.16` and the workspace to `0.0.41` (`Cargo.toml`/`Cargo.lock`).

## v0.0.40

**iso-intron v0.0.15 — `--allow-missing` for tolerant intron lookup**

- Iso-intron v0.0.15 adds `--allow-missing` to handle query introns absent from the reference catalog. When set, both the spliced-intron check (`detect_rt_intron`) and the exon-embedded retention check (`detect_retention`) log a `WARN` and skip the missing `chrom:start-end(strand)` key instead of panicking; without the flag the previous strict behaviour is preserved.
- Fixed the `&Intron` lifetime bug where the `--allow-missing` path tried to return `&Intron` borrowing a temporary synthetic `Intron` (`E0515`). Replaced with skip-and-warn without fabricating a reference, fixing the lifetime and keeping `Schema<'a>` borrowed from the reference map. Threaded the flag from `Args` through `detect_intron_retentions` → `process_components` → `process_component` → `detect_*`.
- Bumped workspace to `0.0.40` and `iso-intron` to `0.0.15` (`Cargo.toml`/`Cargo.lock`).

## v0.0.39

**iso-orphan v0.0.5 — orphan classifier overhaul and per-read TSV report; iso-intron v0.0.14 — `--ignore-utr`; iso-classify v0.0.12 — CDS/UTR intron classification fix**

- Iso-classify v0.0.12 corrects a classification inconsistency for introns that span CDS and UTR regions. An intron already collected in the UTR set is now removed from it the moment it appears in a CDS, enforcing the rule that any intron within a CDS is a CDS intron. Before this fix, such introns could be assigned to the wrong region and distort downstream classification.
- Iso-intron v0.0.14 adds `-u`/`--ignore-utr`, a flag that stops retentions located in UTR regions from being counted as events. When enabled, the affected retentions are excluded from the event tally: a new `W` code (ignored retention) joins the read code legend, the HTML report gains an "Ignored retentions" section, and the run summary now reports how many reads carried ignored introns.
- Iso-orphan v0.0.5 is a major rework of the orphan classification pipeline:
  - Guided and self-guided evaluation are now clearly separated, and in self-guided mode (`--non-overlapping`) the reference step is skipped entirely.
  - Single-exon reads are scored by reciprocal overlap against the best reference exon (`--min-overlap-frac`); multi-exon reads by splice-junction agreement with each reference transcript (`--min-junction-frac`). A rescue is only allowed when the read shares a complete structural feature with one single reference transcript — a lone shared coordinate never rescues a read, and evidence cannot be assembled from mutually incompatible references. Multi-exon reads whose references are all single-exon report `NO_COMPARABLE_REFERENCE_JUNCTIONS` and fall through to de novo evaluation.
  - Component size is treated as context only (`SINGLE_READ_COMPONENT`, `LOW_COMPONENT_SUPPORT`) and never decides an outcome on its own. For multi-exon reads, low abundance now demands stronger evidence instead of causing automatic rejection.
  - De novo evaluation proceeds in three tiers: an intron-chain cluster of at least `--min-read-num-denovo` reads passes unconditionally; failing that, junction support is measured against `--min-intron-support-frac` (a read never counts towards its own junctions); and if structural support is missing or unmeasurable, a median splice-site score reaching `--min-splice-score` can still carry the read. An unavailable splice score never vetoes a read.
  - Single-exon reads keep a terminal abundance requirement through the new `--min-single-exon-support` flag, which defaults to `--min-read-num-denovo`.
  - Splice scoring now compares the **median** score across a read's splice sites instead of the minimum, since BigWig coverage is not uniform and a single uncovered site would otherwise sink an otherwise strong read. The weakest site is still reported as a diagnostic, but it gates nothing.
  - Added `<prefix>.report.tsv`, a per-read report produced on every run: one row per query record explaining which classifier evaluated it, which features and thresholds were compared, what rescued or rejected it, and where its BED record was written. The report is generated from the classification result itself — the same single source of truth that drives BED selection, so the two outputs cannot disagree — and it is never deduplicated, unlike the BED outputs.
  - `--min-discard-percent` is deprecated and has no effect on classification. It is kept only so existing command lines keep parsing, and passing it emits an explicit warning.
  - Fixed a clap constraint that referenced a removed argument, which could panic at startup in debug builds.

## v0.0.38

**iso-intron v0.0.13 — REVIEW suffix preserved alongside classification**

- The `REVIEW` status no longer replaces the existing classification. When a component exceeds the retention ratio threshold during recovery, `/REVIEW` is appended as a suffix (e.g., `PASS/REVIEW`), preserving the original class information so downstream tools can act on both the classification and the review flag independently.
- Added `component_events` and `component_size` fields to the intron component schema, exposing per-component event counts and transcript size directly in the output. These metrics were previously unavailable at the component level.
- HTML report output now includes the component-level breakdown alongside the retention ratio (displayed as `ratio (events/size)`), giving a clearer picture of what drives a component's retention signal.
- Renamed the "RT retentions" section in HTML reports to "RT/Artifact retentions" to accurately reflect what the section covers.

## v0.0.37

**iso-segment v0.0.9 — identity and mapping feature overhaul**

- Rewrote `set_mapping_features` in iso-segment to properly handle bare `M` CIGAR operations (e.g., deSALT output) which caused identity to report 0.0%. The method now correctly counts all aligned columns (`M`, `=`, `X`) together and derives mismatches from the NM tag instead of relying on the CIGAR alone.
- Fixed `end_site` calculation — previously computed as `start + expansion - 1` where *expansion* conflated insertions with reference span. Now it uses `ref_span` (sum of `M`/`=`/`X` + `D` + `N`) to produce correct genomic coordinates.
- Refined clip placement for reverse-strand reads: leading and trailing clips are swapped so that `five_clip` and `three_clip` reflect the biological 5′ and 3′ ends regardless of strand.
- Added a `get_nm` helper to extract the NM tag from BAM records, with safe fallback to `None`.
- Guarded against division by zero and unsigned underflow when alignment length is zero or NM is larger than aligned columns.
- Added a comprehensive test suite covering: M-only CIGAR regression, `=`/`X` equivalence, intron-spanning alignments, indels in the identity denominator, clip placement on both strands, empty alignments, and NM saturation edge cases.

## v0.0.36

**iso-classify and iso-intron — signal strength classes and artifact detection**

- Added `length` branch and `has_minimum_signal` filter to iso-classify for more stringent transcript classification.
- Introduced the `Q` code (`HAS_ARTIFACT`) in iso-intron to flag reads with suspected artifact signals.
- Added three new classification tiers — `HAS_STRONG_RT`, `HAS_WEAK_RT`, and `HAS_ARTIFACT` — providing finer granularity in evaluating retention signal strength.

## v0.0.35

**Missing `--version` flags filled in**

- Added `--version` support to subcommands that were missing it, so every entry point now reports its version when queried.

## v0.0.34

**iso-utr v0.0.10 — duplicate row fix when recovering**

- Fixed a bug in iso-utr where running with `--recover` produced duplicated output rows, ensuring each recovered transcript appears exactly once.

## v0.0.33

**iso-utr v0.0.9 — ghost queries squashed**

- Fixed an issue where iso-utr generated phantom queries during `--recover` operations, which previously inflated output with non-existent records.

## v0.0.32

**iso-align v0.0.4 — back to `--reads` as input, multi-file support**

- Reverted iso-align input from derived FASTA to the original `--reads` flag, which now accepts multiple BAM/SAM files simultaneously for batch alignment.

## v0.0.31

**iso-segment v0.0.8 — fragment tag refinement**

- When `--fragment` is active, iso-segment now writes an additional `F` tag alongside the existing `FG` tag, ensuring unique read names even when fragment re-alignment produces duplicate identifiers.

## v0.0.30

**iso-segment v0.0.7 — fragment re-alignment support**

- Added `-F`/`--fragment` flag to iso-segment, enabling re-alignment of fragment pairs and annotating them with an `FG` tag.
- iso-align v0.0.3 dropped the `--reads` flag in favor of deriving FASTA sequences directly from the input BAM.

## v0.0.29

**iso-segment v0.0.6 — automatic tag correction detection**

- Iso-segment now detects and corrects common tag errors automatically, and the deprecated `--cigar`/`-C` flag has been removed.
- Iso-align v0.0.2 included supplementary alignment handling for more complete read reporting.

## v0.0.28

**iso-align v0.0.1 pre-release**

- First public pre-release of iso-align, a dedicated alignment subcommand for isotools.
- Fixed a bug in iso-fusion's component merging logic.

## v0.0.27

**iso-intron v0.0.11 — whole-chromosome single-exon edge case**

- Patched an edge case where iso-intron would build up an entire chromosome from individual single-exon transcripts, preventing runaway memory and coordinate overflow.

## v0.0.26

**iso-adapter v0.0.3 — BAI default and short flags**

- BAI index generation is now enabled by default.
- Added short option aliases for key adapter-detection flags.

## v0.0.25

**iso-adapter v0.0.2 — poly-A trimming and outward junk detection**

- Implemented `--trim-polya` for trimming poly-A tails from adapter sequences.
- Expanded `--remove-adapters` to catch outward-facing adapter junk, not just inward contamination.

## v0.0.24

**iso-adapter introduced**

- First implementation of iso-adapter, a new subcommand dedicated to adapter detection and removal from long-read transcriptomic data.

## v0.0.23

**iso-segment tag naming fix**

- Changed segment output tags from `{read}{SEP}CG` to `C{name}` format to resolve duplicated read names when processing CIGAR-extended reads.

## v0.0.22

**iso-orphan v0.0.4 and iso-fusion v0.0.12 — junction support and generic reference maps**

- Iso-orphan received a full re-implementation with junction support and splicing score computation.
- Iso-fusion adopted a generic reference map implementation for more flexible fusion detection.

## v0.0.21

**iso-orphan v0.0.3 — scraps and HQ classification**

- Added scraps/high-quality classification to iso-orphan.
- Fixed the output scrap writer and made the reference generic for broader compatibility.

## v0.0.20

**iso-cigar v0.0.2 — alignment extension**

- Added `--extend-alignment` to iso-cigar for extending partial alignments.
- Deleted empty extended alignment entries automatically.
- Enabled BAI index generation by default.

## v0.0.19

**iso-cigar introduced; configuration and packaging overhaul**

- First release of iso-cigar, a CIGAR-based alignment processing subcommand.
- Dropped the `config`, `pack`, `toml` modules and the `iso-split` subcommand to streamline the codebase.
- Added Docker image support and proper documentation infrastructure.

## v0.0.18

**iso-intron NA row handling**

- Allowed iso-intron to process intronIC rows with `NA` values instead of panicking, improving robustness with incomplete annotations.

## v0.0.17

**iso-orphan v0.0.2 introduced**

- First release of iso-orphan, a subcommand for identifying and classifying orphan reads in long-read transcriptomics.

## v0.0.16

**Artifact classification and intron flagging**

- Introduced `Artifact` as a new classification category across iso-intron and iso-classify.
- Added `G` flag in iso-intron v0.0.10 to represent artifact signals.
- Extended classification capabilities to detect and report artifactual features.

## v0.0.15

**iso-fusion v0.0.11 — deduplication and fake output removal**

- Fixed duplicated output rows in iso-fusion fusion reports.
- Removed the separate *fakes* file — fusions previously classified as fake are now sent directly to passes.
- Disabled the `--tag` flag temporarily for stability.

## v0.0.14

**iso-intron v0.0.9 — graceful empty input handling**

- Prevents panic when the input intron set is empty, returning an informative message instead.

## v0.0.13

**Iso-intron gains `--outdir` and prefix adjustments**

- Added `--outdir` to iso-intron for specifying output directories independently.
- Adjusted prefix handling for cleaner output file naming.

## v0.0.12

**iso-utr v0.0.8 re-implemented**

- Full re-implementation of iso-utr with improved UTR detection and annotation logic.

## v0.0.11

**iso-intron stable release; iso-pas introduced; iso-polya dropped**

- Stabilized iso-intron v0.0.8 with a production-ready feature set.
- Introduced iso-pas v0.0.1 for polyadenylation site analysis.
- Dropped the `entry` and `iso-polya` subcommands.

## v0.0.10

**iso-segment v0.0.4 — OOM prevention**

- Implemented streaming dispatch in iso-segment to prevent out-of-memory errors when processing large datasets.

## v0.0.9

**iso-classify v0.0.8 — case-agnostic genome and coordinate fixes**

- Made genome reference loading case-agnostic to handle inconsistent chromosome naming.
- Fixed coordinate extraction from IIC keys used during classification.

## v0.0.8

**iso-classify v0.0.7 — MaxEnt database loading fix**

- Fixed a bug that prevented the MaxEnt splice site scoring database from loading correctly.

## v0.0.7

**iso-classify IIC reading fix**

- Resolved an issue where iso-classify failed to read IIC (isoform-intron classification) input files.

## v0.0.6

**iso-classify v0.0.5 — re-implementation and repeat handling**

- Major re-implementation of iso-classify with improved repeat region handling and classification accuracy.

## v0.0.5

**Pre-release alignment — packbed version bump**

- Bumped packbed version for compatibility with the broader isotools pipeline.

## v0.0.4

**NMD module improvements**

- Updated Cargo lock and dependency versions for the NMD submodule.

## v0.0.3

**iso-fusion v0.0.10 — log cleanup and error expansion**

- Removed misleading log messages from iso-fusion output.
- Expanded error reporting for better debugging of fusion detection.

## v0.0.2

**iso-segment v0.0.2 — split and delimiter support**

- Added `--split` and `--delimiter` options to iso-segment for finer control over segment separation.

## v0.0.1

**Initial release**

- First public release of isotools with core infrastructure.
- Included iso-intron (intron retention analysis), iso-utr (UTR annotation), and packbed (BED file packing utilities).
- Set up CI/CD pipeline, Docker image, and initial project documentation.
