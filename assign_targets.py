#!/usr/bin/env python3
"""Assign MACS3 peaks to target genes through TSS and gene-start windows.
Checked on small artificial test data (example/)."""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

win, lst = lambda x: tuple(int(v) for v in x.split(",")), lambda x: x.split(",")
ap = argparse.ArgumentParser(description=__doc__, epilog="Windows are up,down in bp (negative = upstream); "
                             "pass them as --tss-window=-200,100.")
ap.add_argument("--gff", required=True, help="GFF3 genome annotation")
ap.add_argument("--gff-feature", type=lst, default="gene,CDS", help="GFF feature types to read; a locus with several "
                "rows spans them all")
ap.add_argument("--gff-attribute", default="locus_tag", help="GFF attribute holding the gene id")
ap.add_argument("--peaks", required=True, nargs="+", help="MACS3 narrowPeak file(s), one TF each, named by the file name")
ap.add_argument("--out", required=True, help="output folder; one <TF>_assigned_peaks.bed per peak file, one row per peak "
                "x target gene: chrom, summit start, summit end, target, peak_score, target strand, signal, pscore, qscore, "
                "target_chrom, target_start, target_end, peak_id, operon, anchor (TSS/start_codon), distance_to_anchor")
ap.add_argument("--operons", help="operon BED6; name = operon id or Operon-mapper string")
ap.add_argument("--tss", help="TSS BED6 (TSS = end); name = <tss_id>|<locus_tags 5'->3'>, '.' = unlinked")
ap.add_argument("--no-gene-start", action="store_true", help="no start-codon windows; units without a TSS are not assigned")
ap.add_argument("--no-internal", action="store_true", help="no internal units (operon starts only)")
ap.add_argument("--tss-to-gene", action="store_true", help="a hit through a TSS links only its linked gene")
ap.add_argument("--tss-window", type=win, default="-200,100", help="window around a linked TSS")
ap.add_argument("--start-window", type=win, default="-400,25", help="operon-start window around the start codon")
ap.add_argument("--internal-window", type=win, default="-250,50", help="internal window around the start codon")
a = ap.parse_args()
warn = lambda m: print("warning:", m, file=sys.stderr)

rows, types, attrs = [], set(), set()
for line in open(a.gff):
    f = line.rstrip("\r\n").split("\t")
    if line.startswith("#") or len(f) < 9:
        continue
    types.add(f[2])
    if f[2] in a.gff_feature:
        attrs.update(x.partition("=")[0] for x in f[8].split(";"))
        if m := re.search(rf"(?:^|;){re.escape(a.gff_attribute)}=([^;]+)", f[8]):
            rows.append((m.group(1), f[0], int(f[3]), int(f[4]), f[6]))
g = pd.DataFrame(rows, columns=["locus", "chrom", "start", "end", "strand"])
g = g.groupby("locus", sort=False).agg({"chrom": "first", "start": "min", "end": "max", "strand": "first"}).reset_index()
if g.empty:
    raise SystemExit(f"gff: no {'/'.join(a.gff_feature)} feature with a {a.gff_attribute} attribute. Feature types in the "
                     f"file: {', '.join(sorted(types))}; attributes of the {'/'.join(a.gff_feature)} rows: "
                     f"{', '.join(sorted(attrs)) or 'none'}. Set --gff-feature and --gff-attribute to match.")
sc = dict(zip(g.locus, np.where(g.strand == "+", g.start, g.end)))
strand, chrom = dict(zip(g.locus, g.strand)), dict(zip(g.locus, g.chrom))
g_start, g_end = dict(zip(g.locus, g.start - 1)), dict(zip(g.locus, g.end))  # BED, 0-based start

base = []  # units: genes ordered 5'->3', the first gene's start codon is the unit start
if a.operons:
    op = pd.read_csv(a.operons, sep="\t", header=None, usecols=range(6), dtype={0: str, 3: str},
                     names=["chrom", "start", "end", "name", "score", "strand"])
    gc, gs, ge, gst, gl = (g[k].to_numpy() for k in ("chrom", "start", "end", "strand", "locus"))
    op = op[~op.name.str.contains("type:gene", regex=False)]
    empty = off = 0
    for r in op.itertuples():
        oid = m.group(1) if (m := re.search(r"operon:([^|\]]+)", r.name)) else r.name
        i = np.flatnonzero((gc == r.chrom) & (gst == r.strand) & (gs > r.start) & (ge <= r.end))  # BED start 0-based
        empty += len(i) == 0
        off += len(i) > 0 and not ((gs[i] - 1 == r.start).any() and (ge[i] == r.end).any())
        genes = gl[i[np.argsort(gs[i] * (1 if r.strand == "+" else -1), kind="stable")]].tolist()
        for i in range(min(len(genes), 1 if a.no_internal else len(genes))):
            base.append(dict(operon=oid, internal=bool(i), genes=genes[i:]))
    if empty:
        warn(f"operons: {empty}/{len(op)} rows hold no gene on their strand (chromosome names? 0-based start?)")
    if off:
        warn(f"operons: {off}/{len(op)} rows do not start or end on a gene boundary; edge genes may be lost")
    if lone := (~g.locus.isin({x for u in base for x in u["genes"]})).sum():
        warn(f"operons: {lone}/{len(g)} genes are in no operon and are never assigned")
else:
    base = [dict(operon=".", internal=False, genes=[x]) for x in g.locus]

tss = {}  # linked gene -> [(position, gene list on the TSS's strand and chromosome)]
if a.tss:
    t = pd.read_csv(a.tss, sep="\t", header=None, usecols=range(6), dtype={0: str, 3: str},
                    names=["chrom", "start", "end", "name", "score", "strand"])
    lost = 0
    for r in t.itertuples():
        lst = [x.strip() for x in r.name.partition("|")[2].split(",")]
        same = [x for x in lst if strand.get(x) == r.strand and chrom.get(x) == r.chrom]
        if same and same[0] == lst[0]:
            tss.setdefault(lst[0], []).append((r.end, same))
        elif lst[0] not in ("", "."):
            lost += 1
    if lost:
        warn(f"tss: {lost}/{len(t)} linked rows skipped: gene unknown or on another strand or chromosome")
    if not tss:
        warn("tss: no usable TSS (name must be <tss_id>|<locus_tag>), all units use start-codon windows")

rows = []
for u in base:
    first = u["genes"][0]
    u.update(chrom=chrom[first], strand=strand[first])
    for pos, lst in tss.get(first, []):
        genes = [first] if a.tss_to_gene else u["genes"] if a.operons else lst
        rows.append({**u, "genes": genes, "anchor": "TSS", "anchor_pos": pos, "up": a.tss_window[0],
                     "down": a.tss_window[1], "operon": first if not a.operons and len(genes) > 1 else u["operon"]})
    if first not in tss and not a.no_gene_start:
        lo, hi = a.internal_window if u["internal"] else a.start_window
        rows.append({**u, "anchor": "start_codon", "anchor_pos": sc[first], "up": lo, "down": hi})
units = pd.DataFrame(rows)
if units.empty:
    raise SystemExit("no units to assign")
pos, up, down = units.anchor_pos.to_numpy(), units.up.to_numpy(), units.down.to_numpy()
sgn = np.where(units.strand == "+", 1, -1)
chrom_id = {c: i for i, c in enumerate(units.chrom.unique())}
u_chrom = units.chrom.map(chrom_id).to_numpy()
uid = pd.factorize(units.genes.str[0])[0]  # a unit is identified by its first gene


def resolve(h):
    """Overlapping windows. h: one row per summit p x unit row u whose window holds it (|distance| d, strand s)."""
    h = h.sort_values("d", kind="stable").drop_duplicates(["p", "uid"])  # 1. one anchor per unit: the nearest TSS
    h = h.drop_duplicates(["p", "s"])                       # 2. one unit per strand: the nearest anchor; both strands
    return h.sort_values(["p", "u"])

tfs = [Path(p).stem for p in a.peaks]
if len(set(tfs)) < len(tfs):
    raise SystemExit("peak file names must be distinct: they name the TFs and the output files")
Path(a.out).mkdir(parents=True, exist_ok=True)

for path, tf in zip(a.peaks, tfs):  # one output file per peak file
    try:
        pk = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=range(10), dtype={0: str})
    except pd.errors.EmptyDataError:
        pk = pd.DataFrame({i: [] for i in range(10)}, dtype=int)
    s_all, c_all = (pk[1] + pk[9] + 1).to_numpy(int), pk[0].map(chrom_id).fillna(-1).to_numpy(int)
    if (c_all < 0).any():
        warn(f"{tf}: {(c_all < 0).sum()}/{len(pk)} peaks on chromosomes not in the annotation: "
             f"{sorted(set(pk[0][c_all < 0]))[:3]}")
    out = []
    for c0 in range(0, max(len(pk), 1), 2000):  # chunks bound the size of the site x anchor matrix
        s, c = s_all[c0:c0 + 2000], c_all[c0:c0 + 2000]
        D = (s[:, None] - pos) * sgn
        pi, ui = np.nonzero((D >= up) & (D <= down) & (c[:, None] == u_chrom))
        h = resolve(pd.DataFrame({"p": pi, "u": ui, "d": np.abs(D[pi, ui]), "s": sgn[ui], "uid": uid[ui]}))
        pi, ui = h.p.to_numpy(), h.u.to_numpy()
        u = units.iloc[ui]
        out.append(pd.DataFrame({"i": c0 + pi, "gene": u.genes.to_numpy(), "anchor": u.anchor.to_numpy(),
                                 "distance_to_anchor": D[pi, ui], "operon": u.operon.to_numpy()}))
    r = pd.concat(out, ignore_index=True).explode("gene")
    p = pk.iloc[r.i.to_numpy()]
    bed = pd.DataFrame({
        "chrom": p[0].to_numpy(), "start": s_all[r.i] - 1, "end": s_all[r.i], "target": r.gene.to_numpy(),
        "peak_score": p[4].to_numpy(), "strand": r.gene.map(strand).to_numpy(), "signal": p[6].to_numpy(),
        "pscore": p[7].to_numpy(), "qscore": p[8].to_numpy(), "target_chrom": r.gene.map(chrom).to_numpy(),
        "target_start": r.gene.map(g_start).to_numpy(), "target_end": r.gene.map(g_end).to_numpy(),
        "peak_id": p[3].to_numpy(), "operon": r.operon.to_numpy(), "anchor": r.anchor.to_numpy(),
        "distance_to_anchor": r.distance_to_anchor.to_numpy()})
    bed.sort_values(["chrom", "start"], kind="stable").to_csv(
        Path(a.out) / f"{tf}_assigned_peaks.bed", sep="\t", header=False, index=False)
    print(f"{tf}: {r.i.nunique()}/{len(pk)} peaks assigned, {bed.target.nunique()} target genes")
