# iso-orphan

`iso-orphan` detects orphan reads and background transcription in a query set of long reads. Reads are packed into strand-specific components by genomic overlap, and every query read of a component is classified against two independent sources of evidence.

## Classification

Reference support and independent query support are applied **sequentially**, not as alternative routes:

```text
coherent reference support
    → PASS
otherwise, independent support from the other reads of the component
    → PASS
otherwise
    → SCRAP
```

Reference overlap is therefore evidence, not a routing gate. A read that overlaps a reference exon but fails the reference comparison is not discarded: it continues to de novo evaluation with the same evidence base as a read that had no reference to compare with. A novel isoform at an annotated locus can be retained on the strength of the other reads.

In self-guided mode (`--non-overlapping`) the reference step is skipped entirely.

### Reference support

**Single-exon reads** are scored by reciprocal overlap of their exon with the best reference exon, against `--min-overlap-frac`.

**Multi-exon reads** are scored by splice-junction agreement with each reference transcript individually. If the best query-junction match fraction reaches `--min-junction-frac`, the read is retained. Otherwise a rescue is attempted, and it must be a **complete structural feature shared with one single reference transcript**:

| Evidence | Requirement |
| --- | --- |
| `SPLICE_JUNCTION` | at least one complete junction shared, with every shared junction in the same order in both transcripts |
| `BOTH_TRANSCRIPT_ENDS` | both transcript ends agree with that reference within `--end-tolerance` |
| `EXON_PAIR` | a whole exon agrees with a reference exon at both boundaries |

A single shared coordinate never rescues a read, and a rescue can never be assembled from coordinates contributed by different, mutually incompatible references. Single-exon references are not considered: a spliced read shares no junction with an unspliced transcript, so any agreement is a coincidence rather than support for its structure. A multi-exon read in a component whose references are all single-exon reports `NO_COMPARABLE_REFERENCE_JUNCTIONS` and proceeds to de novo evaluation.

### Independent support

The de novo evidence base is built **once per component, over all of its query reads**. A read's support therefore never depends on how many other reads happened to touch a reference exon.

**Component size is not evidence.** It is a computational grouping property, so it is reported as context (`SINGLE_READ_COMPONENT`, `LOW_COMPONENT_SUPPORT`) but never decides an outcome on its own. For multi-exon reads, low abundance means *stronger evidence is required*, not automatic rejection.

**Multi-exon reads** are evaluated in three tiers:

1. **Intron-chain cluster.** A cluster of at least `--min-read-num-denovo` reads reproducing the same intron chain is a dominant cluster and passes unconditionally.
2. **Junction support.** Otherwise, the fraction of the read's junctions carried by at least `--intron-support-threshold` of the **other** reads is compared against `--min-intron-support-frac`. A read never counts towards the support of its own junctions, and junctions are grouped within `--junction-tolerance`, the same tolerance that clusters intron chains. When support passes, a splice-site model may still veto the read; an *unavailable* score does not veto it.
3. **Splice evidence.** When structural support is missing or unmeasurable, independent sequence evidence can still carry the read: it passes if its median splice-site score reaches `--min-splice-score`. An unavailable score cannot rescue it.

**Single-exon reads** keep a terminal abundance requirement: their reciprocal-overlap cluster must reach `--min-single-exon-support` (default: `--min-read-num-denovo`). They carry no junction, so no orthogonal evidence exists for them.

### Splice-site scores

BigWig scores are 0-based per-base values. For an intron `(start, end)` the two bases carrying splice signal are `start` and `end - 1`; which one is the donor depends on the strand. Each site is read from the track that can carry its signal — donor sites from the donor track, acceptor sites from the acceptor track, on the strand of the read.

Classification compares `--min-splice-score` against the **median** across the read's splice sites. A minimum would be brittle: BigWig coverage is not uniform, so a single uncovered site scores `0.000` and would sink a read whose other sites are all strong. The weakest site is reported alongside as a diagnostic — it shows whether a passing median hides an uncovered or genuinely weak junction — but it gates nothing.

## Outputs

Every run writes three files to `<outdir>/orphans/`:

| File | Content |
| --- | --- |
| `<prefix>.hq.bed` | BED12 records classified as `PASS` |
| `<prefix>.scraps.bed` | BED12 records classified as `SCRAP` |
| `<prefix>.report.tsv` | One row per processed query record, explaining its classification |

The report is always generated; there is no flag to enable it.

## `<prefix>.report.tsv`

### Purpose

The report exposes the classification of every single query read: which classifier evaluated it, which features were calculated, which thresholds were compared, which criteria failed or rescued it, and where its BED record was written.

It is not reconstructed from the BED files. Classification produces one structured result per read, and that result is the single source of truth for the BED selection *and* for the report row, so the two outputs cannot disagree.

Properties worth knowing:

- Exactly one row per valid query BED12 record. Reference records never appear.
- The report is never deduplicated. Query names are not unique, and two byte-identical query records still produce two rows — even though the BED outputs keep their existing deduplication and collapse them into one line.
- Rows are sorted by `chrom, start, end, read_id, decision, reasons`, so the file does not depend on thread scheduling. Running with one thread and with many threads produces the same bytes.

### Columns

| Column | Description |
| --- | --- |
| `read_id` | BED name field, `.` when unavailable |
| `chrom`, `start`, `end`, `strand` | Original BED genomic identity, with BED coordinate semantics (0-based, half-open) |
| `exon_count` | Number of exons of the query record |
| `read_type` | `SINGLE_EXON` or `MULTI_EXON`, derived from the read's intron chain |
| `mode` | Initial program mode: `GUIDED` or `DE_NOVO` |
| `evaluation_path` | Classifier that produced the decision: `GUIDED_REFERENCE`, `GUIDED_TO_DE_NOVO` or `DE_NOVO` |
| `group_size` | Number of query reads in the component, which is the de novo evidence base |
| `reference_count` | Number of reference transcripts in the original packed component |
| `cluster_size` | Size of the intron-chain or reciprocal-overlap cluster assigned to the read, `.` when the read was decided on reference evidence |
| `overlaps_reference_exon` | `true`, `false`, or `.` when there was no reference to test against |
| `best_reference_id` | Reference transcript that produced the best score or carried the rescue, `.` otherwise |
| `best_reference_overlap` | Best reciprocal overlap used for guided single-exon classification, `.` otherwise |
| `junction_matches` | Matched query junctions in the best reference comparison |
| `junction_count` | Total number of query junctions |
| `junction_fraction` | Best query-junction match fraction |
| `boundary_match` | Whether a coherent feature was shared with one reference, `.` when the rescue was not evaluated |
| `boundary_evidence` | Which feature carried the rescue: `SPLICE_JUNCTION`, `BOTH_TRANSCRIPT_ENDS`, `EXON_PAIR`, or `.` |
| `intron_support_fraction` | Fraction of the read's junctions supported by the *other* reads, `.` when not evaluated or unmeasurable |
| `median_splice_score` | Median score across the read's splice sites — the value compared against `--min-splice-score`, `.` when unavailable |
| `weakest_splice_score` | Lowest score across the read's splice sites. Diagnostic only; it gates nothing |
| `applied_thresholds` | Only the thresholds compared on this read's path, `;`-delimited |
| `reasons` | Ordered reason codes, `/`-delimited |
| `decision` | `PASS` or `SCRAP` |

`evaluation_path` is `GUIDED_TO_DE_NOVO` whenever the de novo classifier produced the final decision in guided mode. `overlaps_reference_exon` distinguishes the two ways a read gets there: `false` means it shares no exon with any reference, `true` means its reference evidence was insufficient. A guided component that holds no reference transcript at all reports `evaluation_path=DE_NOVO` with `reference_count=0` and `overlaps_reference_exon=.`, because there is no reference exon to test against.

`intron_support_fraction` is `.` in two distinct situations, told apart by the reason codes: the read passed on reference evidence or on a dominant cluster, so support was never needed, or there was no other read to draw support from (`INTRON_SUPPORT_UNAVAILABLE`). Absent evidence is not evidence of absence.

### Formatting

- `.` marks any unavailable or non-applicable value.
- Booleans are lowercase `true` / `false`.
- Floating-point features and thresholds use three decimal places (`0.500`); integer thresholds are written as integers (`min_cluster_support=5`).
- Reason codes are joined with `/`, applied thresholds with `;`, and neither contains spaces.
- A plain TSV header precedes the data rows; no comment or metadata line is emitted.
- Tabs, carriage returns and newlines inside textual BED fields are escaped (`\t`, `\r`, `\n`), so a single record cannot corrupt the file.

### Decision values

`decision` is exactly `PASS` or `SCRAP`. `PASS` reads are written to `<prefix>.hq.bed`, `SCRAP` reads to `<prefix>.scraps.bed`.

### Reason codes

`reasons` describes the whole evaluation path, not only the branch that produced the decision. A `SCRAP` row lists every applicable failed criterion; a `PASS` row lists the criterion or rescue that justified retention, preceded by the failed criteria that motivated the fallback. Codes appear in evaluation order and are `/`-delimited.

| Code | Meaning |
| --- | --- |
| `REFERENCE_OVERLAP_PASS` | Guided single-exon read met `min_overlap_frac` |
| `LOW_REFERENCE_OVERLAP` | Guided single-exon read below `min_overlap_frac` |
| `JUNCTION_SUPPORT_PASS` | Query-junction match fraction met `min_junction_frac` |
| `LOW_JUNCTION_SUPPORT` | Query-junction match fraction below `min_junction_frac` |
| `NO_COMPARABLE_REFERENCE_JUNCTIONS` | No multi-exon reference was available for junction comparison |
| `BOUNDARY_MATCH_RESCUE` | A complete structural feature shared with one reference rescued the read; see `boundary_evidence` |
| `NO_BOUNDARY_MATCH` | No complete feature was shared with any single reference |
| `NO_REFERENCE_EXON_OVERLAP` | The read shares no exon with any reference |
| `SINGLE_READ_COMPONENT` | The read is alone in its component (context, never terminal) |
| `LOW_COMPONENT_SUPPORT` | The component holds fewer reads than `min_cluster_support` (context, never terminal) |
| `DOMINANT_INTRON_CHAIN_CLUSTER` | The read's intron-chain cluster reaches `min_cluster_support` |
| `LOW_INTRON_CHAIN_CLUSTER_SUPPORT` | The read's intron-chain cluster is below `min_cluster_support` |
| `INTRON_SUPPORT_RESCUE` | The other reads of the component corroborate the read's junctions |
| `LOW_INTRON_SUPPORT` | Junction support from the other reads was below `min_intron_support_frac` |
| `INTRON_SUPPORT_UNAVAILABLE` | No other read carried a junction, so support is unmeasurable |
| `SPLICE_SCORE_PASS` | The median splice-site score met `min_splice_score` after a structural rescue |
| `SPLICE_SCORE_RESCUE` | Splice evidence alone carried a read with no structural support |
| `LOW_SPLICE_SCORE` | The median splice-site score was below `min_splice_score` |
| `SPLICE_SCORE_UNAVAILABLE` | No usable splice score existed: no scores supplied, or the read is unstranded |
| `SUPPORTED_SINGLE_EXON_CLUSTER` | The read's reciprocal-overlap cluster reaches `min_single_exon_support` |
| `LOW_SINGLE_EXON_CLUSTER_SUPPORT` | The read's reciprocal-overlap cluster is below `min_single_exon_support` |

Common combinations:

```text
REFERENCE_OVERLAP_PASS
JUNCTION_SUPPORT_PASS
LOW_JUNCTION_SUPPORT/BOUNDARY_MATCH_RESCUE
LOW_JUNCTION_SUPPORT/NO_BOUNDARY_MATCH/DOMINANT_INTRON_CHAIN_CLUSTER
LOW_REFERENCE_OVERLAP/SUPPORTED_SINGLE_EXON_CLUSTER
NO_REFERENCE_EXON_OVERLAP/DOMINANT_INTRON_CHAIN_CLUSTER
LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_RESCUE/SPLICE_SCORE_PASS
LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_RESCUE/LOW_SPLICE_SCORE
LOW_INTRON_CHAIN_CLUSTER_SUPPORT/LOW_INTRON_SUPPORT/SPLICE_SCORE_RESCUE
SINGLE_READ_COMPONENT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_UNAVAILABLE/SPLICE_SCORE_UNAVAILABLE
```

An unavailable splice score is not treated the same way in the two tiers, and the report states which happened:

- After a **structural rescue**, a missing score does **not** veto the read: `.../INTRON_SUPPORT_RESCUE/SPLICE_SCORE_UNAVAILABLE` is a `PASS`.
- Without structural support, a missing score cannot **rescue** the read: `.../LOW_INTRON_SUPPORT/SPLICE_SCORE_UNAVAILABLE` is a `SCRAP`.

### Example

Reference with exons `[1000,1200)`, `[1400,1600)`, `[1800,2000)`; `read1` reproduces its intron chain, `read2` shares neither a junction nor a whole exon nor both ends with it, and no other read corroborates its junctions:

```text
read_id	chrom	start	end	strand	exon_count	read_type	mode	evaluation_path	group_size	reference_count	cluster_size	overlaps_reference_exon	best_reference_id	best_reference_overlap	junction_matches	junction_count	junction_fraction	boundary_match	boundary_evidence	intron_support_fraction	median_splice_score	weakest_splice_score	applied_thresholds	reasons	decision
read1	chr1	1000	2000	+	3	MULTI_EXON	GUIDED	GUIDED_REFERENCE	2	1	.	true	refA	.	2	2	1.000	.	.	.	.	.	min_junction_frac=0.500	JUNCTION_SUPPORT_PASS	PASS
read2	chr1	1050	1350	+	2	MULTI_EXON	GUIDED	GUIDED_TO_DE_NOVO	2	1	1	true	refA	.	0	1	0.000	false	.	0.000	.	.	min_junction_frac=0.500;junction_tolerance=0;end_tolerance=0;min_cluster_support=5;min_intron_support_frac=0.500;intron_support_threshold=0.500;min_splice_score=0.500	LOW_JUNCTION_SUPPORT/NO_BOUNDARY_MATCH/LOW_COMPONENT_SUPPORT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/LOW_INTRON_SUPPORT/SPLICE_SCORE_UNAVAILABLE	SCRAP
```

`read1` passed on junction support, so no rescue was attempted and no de novo criterion was reached; every column belonging to those steps is `.`. `read2` shows the full sequential path: the reference comparison failed, the coherent rescue found nothing, its component is small, the other reads did not corroborate its junction, and no splice evidence was available to carry it. Note that `LOW_COMPONENT_SUPPORT` appears as context in the middle of the chain — it did not decide anything.

## Usage

```bash
# guided
iso-orphan --query reads.bed --ref annotation.bed --all --outdir results --prefix sample

# de novo
iso-orphan --query reads.bed --non-overlapping --outdir results --prefix sample

# guided, with SpliceAI splice-site scores
iso-orphan --query reads.bed --ref annotation.bed --all --splicing-scores spliceai/ --outdir results
```

### Thresholds

| Option | Applies to |
| --- | --- |
| `--min-junction-frac` | reference junction agreement |
| `--min-overlap-frac` | reference single-exon overlap, and single-exon clustering |
| `--junction-tolerance` | junction matching, intron-chain clustering, junction support grouping |
| `--end-tolerance` | matching both transcript ends to a reference |
| `--min-read-num-denovo` | intron-chain cluster size |
| `--min-single-exon-support` | single-exon overlap cluster size (default: `--min-read-num-denovo`) |
| `--min-intron-support-frac` | fraction of a read's junctions that must be supported |
| `--intron-support-threshold` | fraction of the *other* reads that must carry a junction |
| `--min-splice-score` | median across a read's splice sites |

`--min-discard-percent` is deprecated and has no effect; it was never wired into classification. Run `iso-orphan --help` for the full list.

### Calibrating `--intron-support-threshold`

Junction support is counted leave-one-out, so a read no longer inflates its own support. In a group of *n* reads a junction now needs `threshold × (n - 1)` other reads instead of `threshold × n` reads including itself. Junctions carried by exactly half of a group sit right on the default `0.5` boundary and can fall to the strict side of it, so this threshold is worth re-checking against your data. The `intron_support_fraction` column makes the effect directly measurable.
