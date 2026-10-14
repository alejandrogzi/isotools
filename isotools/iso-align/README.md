# iso-align

Second-pass helper for long introns. minimap2 runs pass 1 with a moderate intron cap
(ARK: `-G 300k`). Reads whose pass-1 alignment shows the shape of an intron that did
not fit are **selected**, realigned with a larger cap (ARK: `-G 1600k`), and the
realignment is **accepted** only when it is better than pass 1.

iso-align does both ends of that loop:

- **selection** (default mode): scan the pass-1 BAM, emit the selected reads as
  FASTA (or names);
- **acceptance** (`--realigned-bam`): compare each realigned read with its pass-1
  alignment and write the names of the reads whose realignment should replace it.

It never edits alignments itself. Replacing records is the caller's job (ARK:
`samtools view -N`, dropping accepted names from the base BAM and keeping only
them in the realigned BAM).

## Selection rules (`--rules`)

Rules are comma-separated; a read is selected when **any** enabled rule fires.
The default is `split` (the original predicate). ARK runs all four.

Selection must read the **pass-1 BAM**, before iso-cigar: iso-cigar stitches split
primaries across the gap and drops their `SA` tag, which hides them from `split`.

| rule | records read | fires when | what it catches |
|---|---|---|---|
| `split` | all primary + supplementary records of the read | ≥ 2 records, same chromosome and strand; the largest gap between consecutive records (sorted by start) is in `[--min-gap, --max-gap]`; and a clip ≥ `--min-flank-clip` sits on the gap-facing flank (`--flank-side inner`) or on any end (`any`) | an exon beyond the intron cap, placed as a supplementary alignment |
| `insertion` | primary only | an insertion ≥ `--ins-min-len` (20 bp) lies within `--ins-max-dist` (10 bp) of an `N` ≥ `--ins-min-skip` (20 kb). The distance is from the insertion to the nearer boundary of its closest `N` | a **bridge**: an exon left unaligned (as an insertion) between two long introns, so minimap2 bridged the gene with one intron |
| `long-intron-clip` | primary only | the longest `N` ≥ `--long-intron-min-skip` (150 kb) and a clip ≥ `--long-intron-min-clip` (21 bp, i.e. > 20) on either end | a terminal exon past a long intron, clipped off instead of split off |
| `orphan-clip` | primary only | no `SA` tag and a clip ≥ `--orphan-min-clip` (50 bp) on either end | a terminal exon minimap2 could not place at all |

Clips count soft and hard clips. Secondary and unmapped records are ignored.

ARK settings: `--min-gap 150000 --max-gap 1600000 --min-flank-clip 20 --flank-side inner
--rules split,insertion,long-intron-clip,orphan-clip`.

## Acceptance (`--realigned-bam`)

`--bam` must be the **pass-1** BAM, because iso-cigar drops `AS`/`NM` on the records it
edits. Only primary records are compared. A selected read is **accepted** when all of
the following hold, checked in this order:

1. **Scores present:** both primaries carry `AS` (otherwise `missing_score`, rejected).
2. **Better score:** `AS(realigned) > AS(pass-1)`, strictly (otherwise `rejected_score`).
3. **High quality:** the realigned primary has identity ≥ `--min-identity` (98.0) and a
   clip ≤ `--max-clip` (20 bp) on **both** ends (otherwise `rejected_quality`). A
   realigned primary without `NM` has no identity and is rejected.

Identity is computed as iso-segment computes it, so an accepted read is also an
iso-segment hq read:

```
aligned    = M + = + X
indels     = I + D
mismatches = NM - indels
identity   = round(1000 * (aligned - mismatches) / (aligned + indels)) / 10
```

Output: the sorted names of accepted reads, one per line. `--report` writes, per read,
the pass-1 / realigned `AS`, identity, clips and the decision. A missing identity
prints as `NA`.

### Why the quality gate

Without it, every selected read was replaced blindly. The sweep over identity
thresholds 90 to 99 on the simulation, compared with blind replacement of the
`split` reads only (harmed = exact chain with blind `split`, not exact after;
fixed = the reverse), gave:

| threshold | harmed | fixed |
|---|---|---|
| ≤ 96 | 0 | 227 |
| 97 | 11 | 226 |
| 98 | 51 | 212 |
| 99 | 250 | 187 |

On real data (10 mice), known long-intron recall was flat across thresholds, and 98
added the fewest unsupported long introns. 98 is the default: it matches
iso-segment's hq cut, so accepted reads are not dropped later.

## Examples

```sh
# selection: FASTA of reads to realign
iso-align --bam pass1.bam --reads chunk.fasta.gz -o fragments.fasta \
  --min-gap 150000 --max-gap 1600000 --rules split,insertion,long-intron-clip,orphan-clip

# acceptance: names of reads whose realignment replaces pass 1
iso-align --bam pass1.bam --realigned-bam fragments.bam -o accepted.txt \
  --min-identity 98 --max-clip 20 --report accept.tsv
```
