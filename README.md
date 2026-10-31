# assign_targets.py

Assigns transcription factor (TF) binding peaks from ChIP-seq or DAP-seq, called with MACS3, to the genes they
probably regulate. A peak is assigned to a gene when its summit (the single base with the highest signal) lies in
a window around that gene's transcription start site (TSS), or, where no TSS is known, around its start codon.
Genes transcribed together in one operon are assigned together.

Status: the default mode (operons + TSS map + start-codon fallback), runs
without an operon file or without a TSS map, the options `--no-gene-start`, `--no-internal` and `--tss-to-gene`,
plain operon ids and genomes with two chromosomes give the expected results
on small artificial test data.

## Requirements

Python ≥ 3.8 with numpy and pandas.

## Usage

```bash
python assign_targets.py --gff genome.gff --peaks peaks/*.narrowPeak --operons operons.bed --tss tss.bed --out out/
```

Only `--gff`, `--peaks` and `--out` are required. Each peak file is processed in turn and gives one output file.

## Example

`example/` holds an artificial genome to show the file formats: 100 genes on two chromosomes (`chrom1`,
`plasmid1`), 45 operons, 28 TSSs and peak files for four made-up TFs. `example/expected/` is the output of

```bash
cd example
python ../assign_targets.py --gff genome.gff --operons operons.bed --tss tss.bed --peaks peaks/*.narrowPeak --out out
```

which should match `expected/` file for file.

## Input files

All coordinates refer to the same genome; the chromosome name (column 1, or the GFF seqid) must be written
identically in every file. Genomes with several chromosomes or plasmids are supported.
Rows that cannot be matched to the annotation (peaks on unknown chromosomes, operons holding no gene or not
ending on a gene boundary, linked TSSs whose gene is not found, genes in no operon) are counted in warnings on stderr.

### `--gff` genome annotation (required)

GFF3. Rows of the listed feature types are read (`--gff-feature`, default `gene,CDS`); a locus with several rows,
of one type or of several, spans all of them. Every row needs the attribute holding the gene identifier
(`--gff-attribute`, default `locus_tag`), which is used in all other files and in the output.

| column | used |
|---|---|
| 1 seqid | chromosome |
| 3 type | one of `--gff-feature` (default `gene,CDS`) |
| 4 start, 5 end | 1-based, inclusive |
| 7 strand | `+` or `-`; the start codon is `start` on `+` and `end` on `-` |
| 9 attributes | `<--gff-attribute>=...` |

### `--peaks` peak files (required, one or more)

MACS3 narrowPeak: 10 tab-separated columns, no header, one file per TF.

| column | content |
|---|---|
| 1 chrom | chromosome |
| 2 start | 0-based peak start |
| 3 end | peak end |
| 4 name | peak id (copied to the output) |
| 5 score | copied to the output |
| 6 strand | ignored |
| 7 signal | fold enrichment, copied to the output |
| 8 pValue, 9 qValue | −log10, copied to the output |
| 10 peak | summit offset from column 2 |

The summit is column 2 + column 10 (0-based). The TF name is the file name without extension, so file names must
be distinct. An empty peak file gives an empty output file.

### `--operons` operon map (optional)

BED6, tab-separated, no header, one row per operon.

| column | content |
|---|---|
| 1 chrom | chromosome |
| 2 start | 0-based start of the operon's leftmost gene (GFF start − 1) |
| 3 end | end of the operon's rightmost gene (= GFF end) |
| 4 name | operon id, or an Operon-mapper string `[type:operon\|operon:<id>\|...]` (id = the `operon:` field); rows containing `type:gene` are skipped |
| 5 score | ignored |
| 6 strand | `+` or `-` |

The operon's genes are the GFF genes on the same chromosome and strand that lie entirely between start and end.

### `--tss` TSS map (optional)

BED6, tab-separated, no header, one row per TSS.

| column | content |
|---|---|
| 1 chrom | chromosome |
| 2 start, 3 end | the TSS as a 1-bp interval (0-based start; TSS position = end, 1-based) |
| 4 name | `<tss_id>\|<locus_tags>`: comma-separated locus tags, 5'→3', of the genes this TSS transcribes; `.` = not linked to a gene |
| 5 score | ignored |
| 6 strand | `+` or `-` |

The first locus tag is the linked gene, the gene the TSS starts. The others are the genes after it in the same
transcript, if known; they are used only without `--operons`. TSSs that are unlinked, or whose linked gene lies
on another strand or chromosome, are skipped, as are listed genes on another strand or chromosome.

Example: `NC_000913.3  336  337  TSS0001|b0001,b0002,b0003  0  +`

## Assignment

**Units.** With `--operons`, every operon is a unit (operon start), and every gene of an operon other than the
first starts an internal unit: that gene and the genes downstream of it. Genes in no operon of the file are never
assigned, so the file must list single-gene operons too. Without `--operons`, every gene is a
unit; a TSS hit assigns the genes of its TSS list.

**Windows** (bp, along the unit's strand, negative = upstream; summit position minus anchor position):

| unit | anchor | default | option |
|---|---|---|---|
| any unit with a linked TSS | the TSS | −200 to +100 | `--tss-window` |
| operon start without a TSS | start codon of the first gene | −400 to +25 | `--start-window` |
| internal start without a TSS | start codon of that gene | −250 to +50 | `--internal-window` |

Pass windows with `=`, as in `--tss-window=-150,100`, so the leading minus is not read as an option.

**Overlapping windows** (function `resolve` in the script).
1. A unit with several TSSs is assigned once, through the TSS nearest the summit.
2. Of several units on one strand, only the one with the nearest anchor is assigned.
3. A summit in windows on both strands (between divergent genes) is assigned on both.
4. An assigned unit assigns all its genes.

## Options

| option | effect |
|---|---|
| `--gff-feature` | GFF feature types to read, comma-separated (default `gene,CDS`); a locus spans all its rows |
| `--gff-attribute` | GFF attribute holding the gene id (default `locus_tag`) |
| `--no-gene-start` | no start-codon windows; units without a TSS are not assigned |
| `--no-internal` | operon starts only, no internal units |
| `--tss-to-gene` | a hit through a TSS assigns only the gene the TSS starts |
| `--tss-window`, `--start-window`, `--internal-window` | window `up,down` in bp |

## Output

One `<TF>_assigned_peaks.bed` per peak file in `--out`: tab-separated, no header, one row per peak and assigned
gene, sorted by position. A peak assigned to an operon of four genes gives four rows.

| column | content |
|---|---|
| 1 chrom | chromosome of the peak |
| 2 start, 3 end | summit as a 1-bp interval (0-based start) |
| 4 target | locus tag of the assigned gene |
| 5 peak_score | narrowPeak column 5 |
| 6 strand | strand of the target gene |
| 7 signal, 8 pscore, 9 qscore | narrowPeak columns 7–9 |
| 10 target_chrom | chromosome of the target gene |
| 11 target_start, 12 target_end | target gene coordinates (0-based start) |
| 13 peak_id | narrowPeak column 4 |
| 14 operon | operon id; without `--operons`, the TSS's linked gene when its list holds several genes, else `.` |
| 15 anchor | `TSS` or `start_codon` |
| 16 distance_to_anchor | summit minus anchor, along the target's strand (negative = upstream) |

The script prints one line per TF: peaks assigned and number of target genes.
