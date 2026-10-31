# assign_targets.py

Assigns transcription factor (TF) binding peaks from ChIP-seq or DAP-seq, called with MACS3, to the genes they
probably regulate. A peak is assigned to a gene when its summit (the single base with the highest signal) lies in
a window around that gene's transcription start site (TSS), or, where no TSS is known, around its start codon.
Genes transcribed together in one operon are assigned together. Without an operon map, `predict_operons.py`
predicts operons from gene spacing.

Status: the default mode (operons + TSS map + start-codon fallback), runs
without an operon file or without a TSS map, the options `--no-gene-start` and `--no-internal`, and genomes with
two chromosomes give the expected results on small artificial test data.

## Requirements

Python ≥ 3.8 with numpy and pandas.

## Usage

```bash
python assign_targets.py --gff genome.gff --peaks peaks/*.narrowPeak --operons operons.tsv --tss tss.bed --out out/
```

Only `--gff`, `--peaks` and `--out` are required. Each peak file is processed in turn and gives one output file.

The script can also be imported: `read_genome`, `build_units`, `read_peaks` and `assign` take paths or open files and
return tables and warning lists; `main()` is the command line.

## Example

`example/` holds an artificial genome to show the file formats: 100 genes on two chromosomes (`chrom1`,
`plasmid1`), 45 operons, 28 TSSs and peak files for four made-up TFs. `example/expected/` is the output of

```bash
cd example
python ../assign_targets.py --gff genome.gff --operons operons.tsv --tss tss.bed --peaks peaks/*.narrowPeak --out out
```

which should match `expected/` file for file.

## Input files

All coordinates refer to the same genome; the chromosome name (column 1, or the GFF seqid) must be written
identically in every file. Genomes with several chromosomes or plasmids are supported. Files are plain text, named
by format: `.gff` or `.gff3`, `.tsv`, `.bed`, `.narrowPeak`; a file of the wrong shape (columns, header lines, non-integer
coordinates, missing summits) stops the run with a message naming it.
Rows that cannot be matched to the annotation (loci with rows on several chromosomes or strands, peaks on unknown
chromosomes, operon genes not in the GFF, operons with genes on several chromosomes or strands, linked
TSSs whose gene is not found) are counted in warnings on stderr.

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

One file per TF, either MACS3 narrowPeak (`.narrowPeak`, 10 tab-separated columns, no header) or a BED file of
summits (`.bed`, every row a 1-bp interval at the summit, as in MACS3 `NAME_summits.bed`; columns 4 and 5 are
copied to the output, signal, pscore and qscore are written as `.`).

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

The summit is column 2 + column 10 (0-based), or column 3 of a summit BED. The TF name is the file name without
extension, so file names must
be distinct. An empty peak file gives an empty output file.

### `--operons` operon map (optional)

Tab-separated, no header, two columns: operon id, locus tags of its genes separated by commas. One row per operon;
rows with the same operon id are merged, so one gene per row works too. The order of the genes in the list does
not matter: the genes of an operon are ordered 5'→3' from their GFF coordinates. Genes not listed are single-gene
units. Listed genes that are not in the GFF, and operons whose genes lie on several chromosomes or strands, are
skipped with a warning.

```
op1	b0002,b0003,b0004
op2	b0005
```

### `--tss` TSS map (optional)

BED6, tab-separated, no header, one row per TSS.

| column | content |
|---|---|
| 1 chrom | chromosome |
| 2 start, 3 end | the TSS as a 1-bp interval (0-based start; TSS position = end, 1-based) |
| 4 name | locus tag of the gene the TSS starts; `.` = not linked to a gene |
| 5 score | ignored |
| 6 strand | `+` or `-` |

A TSS is an anchor for the unit its gene starts: the operon or internal unit beginning with that gene, or the gene
itself. TSSs that are unlinked, or whose gene is not in the GFF or lies on another strand or chromosome, are skipped.

Example: `NC_000913.3  336  337  b0001  0  +`

## Assignment

**Units.** With `--operons`, every operon is a unit (operon start), and every gene of an operon other than the
first starts an internal unit: that gene and the genes downstream of it. Genes in no operon are single-gene
units. Without `--operons`, every gene is a unit.

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

## Predicted operons

Without an operon map, `predict_operons.py` can predict operons from gene spacing using the method from UniOP (Su,
Zhang and Söding, bioRxiv 2024, doi:10.1101/2024.11.11.623000). Pairs at or above the cutoff (default 0.5) are joined;
operons of two or more genes are written in the operon map format.

```bash
python predict_operons.py --gff genome.gff --out operons.tsv [--pairs pairs.tsv] [--cutoff 0.5]
```

`--gff-feature` and `--gff-attribute` work as in `assign_targets.py`. The pair table has one row per pair of
consecutive genes:

| column | content |
|---|---|
| chrom | chromosome |
| a, b | locus tags, in genome order |
| gap | start of b − end of a + 1 |
| class | `same`, `convergent` or `divergent` |
| prob | probability that a and b are in one operon; empty for opposite-strand pairs |

The probabilities equal those of the UniOP script (`tests/`), except that pairs are formed within a chromosome
only, the prior is clipped to 0.01–0.99, and a genome without convergent or without divergent pairs uses the
other class alone. The method needs a genome-sized gene set: below about 1,000 same-strand pairs the gaps stop
separating the pairs and nearly all same-strand neighbours are joined.

## Options

| option | effect |
|---|---|
| `--gff-feature` | GFF feature types to read, comma-separated (default `gene,CDS`); a locus spans all its rows |
| `--gff-attribute` | GFF attribute holding the gene id (default `locus_tag`) |
| `--no-gene-start` | no start-codon windows; units without a TSS are not assigned |
| `--no-internal` | operon starts only, no internal units |
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
| 14 operon | operon id, `.` for a gene in no operon |
| 15 anchor | `TSS` or `start_codon` |
| 16 distance_to_anchor | summit minus anchor, along the target's strand (negative = upstream) |

The script prints one line per TF: peaks assigned and number of target genes.
