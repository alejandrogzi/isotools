# iso-cigar

`iso-cigar` rescues missed 3' splice junctions from BAM alignments by converting terminal 3' soft-clipped sequence into downstream exon matches when the sequence agrees with the reference and the annotation.

## Inputs

- A BAM file.
- An annotation file (`BED12`, `GTF`, or `GFF`) used to enumerate transcript exon chains.
- A reference genome in `2bit`.

## Annotation Index

At load time the tool builds two linked annotation views:

- A deduplicated boundary index keyed by `(contig, strand, boundary)` for fast candidate lookup.
- Transcript-aware models that preserve exon order and the mapping from each starting junction to the transcripts that support it.

## Read Selection

For each primary alignment:

1. Reorient the CIGAR and query into transcript 5'->3' order.
2. Require a terminal 3' soft clip with length `>= --clip-cutoff`.
3. Extract a short matched context immediately upstream of the clip.
4. Query annotation boundaries within `+- --wiggle`.

Reads that fail any of these checks are emitted unchanged.

## Rescue

For each candidate starting junction:

1. Compute the strand-normalized boundary delta.
2. If the read stops upstream of the annotated boundary, validate the missing upstream suffix against the upstream exon.
3. Match the clipped query against the first downstream exon.
4. Without `--extend-alignment`, stop here.
5. With `--extend-alignment`, walk downstream exon-by-exon only along transcript-supported paths and keep extending while the query matches exactly.

The walk stops on query exhaustion, exon exhaustion, immediate mismatch, or a partial exon match.

## Rewrite

If rescue is valid:

- Adjust the terminal aligned block by `delta`.
- Replace rescued soft clip bases with one or more `N` + `=` blocks in transcript order.
- Reverse the op list back for reverse-strand genomic BAM output.
- Recompute reverse-strand POS when needed.
- Remove stale alignment tags and attach rescue metadata tags.

## Candidate Choice

Distinct corrected alignments are deduplicated by `(start, cigar)` and scored by:

1. rescued bases taken from the soft clip,
2. number of downstream exons rescued,
3. smaller absolute `delta`,
4. deterministic path order.

By default only the best corrected primary is emitted. With `--keep-additional-corrections`, geometrically distinct lower-ranked corrections are emitted as secondary alignments.

## Outputs

- Without `--split-bam`, all output goes to `<input>.extended.bam`.
- With `--split-bam`, corrected records go to `<input>.extended.bam` and unchanged records go to `<input>.aligned.bam`.
- A `.bai` is written for each emitted BAM as `<bam>.bai`.
- If `--split-bam` is set and no primary record was corrected, the empty `.extended.bam` and `.extended.bam.bai` are removed at the end.
