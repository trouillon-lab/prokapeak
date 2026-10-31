#!/usr/bin/env python3
"""Assign MACS3 peaks to target genes through TSS and gene-start windows.
Checked on small artificial test data (example/)."""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from predict_operons import CUTOFF, predict_operons

TSS_WINDOW, START_WINDOW, INTERNAL_WINDOW = (-200, 100), (-400, 25), (-250, 50)  # up, down in bp around the anchor
CHUNK = 50  # peaks per pass; bounds the size of the summit x anchor matrix


def ext(path, *ok):
    if Path(path).suffix not in ok:
        raise ValueError(f"{path}: file name must end in {' or '.join(ok)} (plain text, not gzipped)")
    return Path(path).suffix


def read_table(path, n, what, num=(1, 2, 9)):
    """Tab-separated, no header, at least n columns; columns num (0-based, those < n) must be integers."""
    try:
        t = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    except pd.errors.EmptyDataError:
        t = pd.DataFrame(columns=range(n), dtype=str)
    except pd.errors.ParserError:
        raise ValueError(f"{path}: {what}, tab-separated without header lines")
    if t.shape[1] < n:
        raise ValueError(f"{path}: {what}; found {t.shape[1]} columns")
    num = t[[c for c in num if c < n]].apply(pd.to_numeric, errors="coerce")
    if num.isna().any().any():
        raise ValueError(f"{path}: {what}; non-integer coordinates (header line?)")
    t[num.columns] = num.astype(int)
    return t


def read_genome(gff, features=("gene", "CDS"), attribute="locus_tag"):
    """GFF3 path or file -> one span per locus (locus, chrom, start, end, strand), list of warnings."""
    rows, types, attrs, w = [], set(), set(), []
    for line in (gff if hasattr(gff, "read") else open(gff)):
        f = line.rstrip("\r\n").split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        types.add(f[2])
        if f[2] in features:
            attrs.update(x.partition("=")[0] for x in f[8].split(";"))
            if m := re.search(rf"(?:^|;){re.escape(attribute)}=([^;]+)", f[8]):
                rows.append((m.group(1), f[0], int(f[3]), int(f[4]), f[6]))
    g = pd.DataFrame(rows, columns=["locus", "chrom", "start", "end", "strand"])
    split = g.groupby("locus")[["chrom", "strand"]].nunique().max(axis=1) > 1
    if split.any():
        w.append(f"gff: {split.sum()} loci with rows on several chromosomes or strands dropped: {list(split.index[split][:3])}")
        g = g[~g.locus.isin(split.index[split])]
    g = g.groupby("locus", sort=False).agg({"chrom": "first", "start": "min", "end": "max", "strand": "first"}).reset_index()
    if g.empty:
        raise ValueError(f"gff: no {'/'.join(features)} feature with a {attribute} attribute. Feature types in the "
                         f"file: {', '.join(sorted(types))}; attributes of the {'/'.join(features)} rows: "
                         f"{', '.join(sorted(attrs)) or 'none'}. Set --gff-feature and --gff-attribute to match.")
    return g, w


def build_units(g, operons=None, tss=None, tss_window=TSS_WINDOW, start_window=START_WINDOW,
                internal_window=INTERNAL_WINDOW, no_internal=False, no_gene_start=False):
    """Units to assign (operon starts, internal starts or single genes) with anchor and window, list of warnings.
    operons: table with operon id and comma-separated locus tags per row; tss: BED6 table; from read_table, or None."""
    w = []
    sc = dict(zip(g.locus, np.where(g.strand == "+", g.start, g.end)))
    strand, chrom = dict(zip(g.locus, g.strand)), dict(zip(g.locus, g.chrom))
    base = []  # units: genes ordered 5'->3', the first gene's start codon is the unit start
    if operons is not None:
        op = operons.iloc[:, :2].set_axis(["operon", "locus"], axis=1)
        op = op.assign(locus=op.locus.str.strip().str.split(r"\s*,\s*", regex=True)).explode("locus")  # one row per gene
        if unknown := (~op.locus.isin(g.locus)).sum():
            w.append(f"operons: {unknown}/{len(op)} listed genes are not in the annotation")
        order = dict(zip(g.locus, np.where(g.strand == "+", g.start, -g.start)))  # 5'->3' along the strand
        split = 0
        for oid, genes in op[op.locus.isin(g.locus)].groupby("operon", sort=False).locus:
            genes = sorted(genes, key=order.get)
            if len({(chrom[x], strand[x]) for x in genes}) > 1:
                split += 1
                continue
            for i in range(1 if no_internal else len(genes)):
                base.append(dict(operon=oid, internal=bool(i), genes=genes[i:]))
        if split:
            w.append(f"operons: {split} operons with genes on several chromosomes or strands ignored")
    listed = {x for u in base for x in u["genes"]}
    base += [dict(operon=".", internal=False, genes=[x]) for x in g.locus if x not in listed]  # single-gene units

    links = {}  # gene -> positions of the TSSs that start it, on its strand and chromosome
    if tss is not None:
        t = tss.iloc[:, :6].set_axis(["chrom", "start", "end", "name", "score", "strand"], axis=1)
        ok = (t["name"].map(strand) == t.strand) & (t["name"].map(chrom) == t.chrom)
        for r in t[ok].itertuples():
            links.setdefault(r.name, []).append(r.end)
        if lost := (~ok & ~t["name"].isin(["", "."])).sum():
            w.append(f"tss: {lost}/{len(t)} linked rows skipped: gene unknown or on another strand or chromosome")
        if not links:
            w.append("tss: no usable TSS (name must be the locus tag of the gene it starts), all units use start-codon windows")

    rows = []
    for u in base:
        first = u["genes"][0]
        u.update(chrom=chrom[first], strand=strand[first])
        for pos in links.get(first, []):
            rows.append({**u, "anchor": "TSS", "anchor_pos": pos, "up": tss_window[0], "down": tss_window[1]})
        if first not in links and not no_gene_start:
            lo, hi = internal_window if u["internal"] else start_window
            rows.append({**u, "anchor": "start_codon", "anchor_pos": sc[first], "up": lo, "down": hi})
    return pd.DataFrame(rows), w


def resolve(h):
    """Overlapping windows: every rule deciding which window a summit is assigned to lives here.
    h: one row per summit p x unit row u whose window holds it (|distance| d, strand s)."""
    h = h.sort_values("d", kind="stable").drop_duplicates(["p", "uid"])  # 1. one anchor per unit: the nearest TSS
    h = h.drop_duplicates(["p", "s"])                       # 2. one unit per strand: the nearest anchor; both strands
    return h.sort_values(["p", "u"])


def read_peaks(path, fmt):
    """narrowPeak (.narrowPeak) or 1-bp summit BED (.bed) -> 10-column table, 1-based summit positions."""
    if fmt == ".narrowPeak":
        pk = read_table(path, 10, "narrowPeak: 10 columns, summit offset in column 10")
        if (pk[9] < 0).any():
            raise ValueError(f"{path}: {(pk[9] < 0).sum()} rows without summit (column 10 = -1)")
        return pk, (pk[1] + pk[9] + 1).to_numpy(int)
    pk = read_table(path, 3, "summit BED: 1-bp intervals at the summits")
    if (pk[2] - pk[1] != 1).any():
        raise ValueError(f"{path}: {(pk[2] - pk[1] != 1).sum()} rows are not 1-bp summit intervals")
    return pk.reindex(columns=range(10), fill_value="."), pk[2].to_numpy(int)


def assign(units, g, pk, summits):
    """One TF: peak table and summits -> output table (one row per peak x target gene), peaks assigned, warnings."""
    if units.empty:
        raise ValueError("no units to assign")
    pos, up, down = units.anchor_pos.to_numpy(), units.up.to_numpy(), units.down.to_numpy()
    sgn = np.where(units.strand == "+", 1, -1)
    chrom_id = {c: i for i, c in enumerate(units.chrom.unique())}
    u_chrom = units.chrom.map(chrom_id).to_numpy()
    uid = pd.factorize(units.genes.str[0])[0]  # a unit is identified by its first gene
    c_all, w = pk[0].map(chrom_id).fillna(-1).to_numpy(int), []
    if (c_all < 0).any():
        w.append(f"{(c_all < 0).sum()}/{len(pk)} peaks on chromosomes not in the annotation: "
                 f"{sorted(set(pk[0][c_all < 0]))[:3]}")
    out = []
    for c0 in range(0, max(len(pk), 1), CHUNK):
        s, c = summits[c0:c0 + CHUNK], c_all[c0:c0 + CHUNK]
        D = (s[:, None] - pos) * sgn
        pi, ui = np.nonzero((D >= up) & (D <= down) & (c[:, None] == u_chrom))
        h = resolve(pd.DataFrame({"p": pi, "u": ui, "d": np.abs(D[pi, ui]), "s": sgn[ui], "uid": uid[ui]}))
        pi, ui = h.p.to_numpy(), h.u.to_numpy()
        u = units.iloc[ui]
        out.append(pd.DataFrame({"i": c0 + pi, "gene": u.genes.to_numpy(), "anchor": u.anchor.to_numpy(),
                                 "distance_to_anchor": D[pi, ui], "operon": u.operon.to_numpy()}))
    r = pd.concat(out, ignore_index=True).explode("gene")
    p = pk.iloc[r.i.to_numpy()]
    strand, chrom = dict(zip(g.locus, g.strand)), dict(zip(g.locus, g.chrom))
    g_start, g_end = dict(zip(g.locus, g.start - 1)), dict(zip(g.locus, g.end))  # BED, 0-based start
    bed = pd.DataFrame({
        "chrom": p[0].to_numpy(), "start": summits[r.i] - 1, "end": summits[r.i], "target": r.gene.to_numpy(),
        "peak_score": p[4].to_numpy(), "strand": r.gene.map(strand).to_numpy(), "signal": p[6].to_numpy(),
        "pscore": p[7].to_numpy(), "qscore": p[8].to_numpy(), "target_chrom": r.gene.map(chrom).to_numpy(),
        "target_start": r.gene.map(g_start).to_numpy(), "target_end": r.gene.map(g_end).to_numpy(),
        "peak_id": p[3].to_numpy(), "operon": r.operon.to_numpy(), "anchor": r.anchor.to_numpy(),
        "distance_to_anchor": r.distance_to_anchor.to_numpy()})
    return bed.sort_values(["chrom", "start"], kind="stable"), r.i.nunique(), w


def main(argv=None):
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
    ap.add_argument("--operons", help="operon map: tab-separated, operon id and its genes (comma-separated), one row per operon")
    ap.add_argument("--tss", help="TSS BED6 (TSS = end); name = locus tag of the gene the TSS starts, '.' = unlinked")
    ap.add_argument("--predict-operons", action="store_true", help="predict operons from gene spacing instead of --operons; "
                    "writes predicted_operons.tsv and predicted_operon_pairs.tsv into --out")
    ap.add_argument("--operon-cutoff", type=float, default=CUTOFF, help="pair probability that joins two genes (--predict-operons)")
    ap.add_argument("--no-gene-start", action="store_true", help="no start-codon windows; units without a TSS are not assigned")
    ap.add_argument("--no-internal", action="store_true", help="no internal units (operon starts only)")
    ap.add_argument("--tss-window", type=win, default="%d,%d" % TSS_WINDOW, help="window around a linked TSS")
    ap.add_argument("--start-window", type=win, default="%d,%d" % START_WINDOW, help="operon-start window around the start codon")
    ap.add_argument("--internal-window", type=win, default="%d,%d" % INTERNAL_WINDOW, help="internal window around the start codon")
    a = ap.parse_args(argv)
    warn = lambda m: print("warning:", m, file=sys.stderr)
    try:
        ext(a.gff, ".gff", ".gff3")
        g, w = read_genome(a.gff, a.gff_feature, a.gff_attribute)
        if a.predict_operons and a.operons:
            raise ValueError("--predict-operons cannot be combined with --operons")
        op = read_table(a.operons, 2, "operon file: operon id, locus tags", num=()) if a.operons and ext(a.operons, ".tsv") else None
        if a.predict_operons:
            op, pairs, w3 = predict_operons(g, a.operon_cutoff)
            Path(a.out).mkdir(parents=True, exist_ok=True)
            op.to_csv(Path(a.out) / "predicted_operons.tsv", sep="\t", header=False, index=False)
            pairs.to_csv(Path(a.out) / "predicted_operon_pairs.tsv", sep="\t", index=False, float_format="%.6f")
            w += w3
        ts = read_table(a.tss, 6, "TSS file: BED6") if a.tss and ext(a.tss, ".bed") else None
        units, w2 = build_units(g, op, ts, a.tss_window, a.start_window, a.internal_window, a.no_internal, a.no_gene_start)
        for m in w + w2:
            warn(m)
        if units.empty:
            raise ValueError("no units to assign")
        tfs, fmts = [Path(p).stem for p in a.peaks], [ext(p, ".narrowPeak", ".bed") for p in a.peaks]
        if len(set(tfs)) < len(tfs):
            raise ValueError("peak file names must be distinct: they name the TFs and the output files")
        Path(a.out).mkdir(parents=True, exist_ok=True)
        for path, tf, fmt in zip(a.peaks, tfs, fmts):  # one output file per peak file
            pk, summits = read_peaks(path, fmt)
            bed, n, w3 = assign(units, g, pk, summits)
            for m in w3:
                warn(f"{tf}: {m}")
            bed.to_csv(Path(a.out) / f"{tf}_assigned_peaks.bed", sep="\t", header=False, index=False)
            print(f"{tf}: {n}/{len(pk)} peaks assigned, {bed.target.nunique()} target genes")
    except ValueError as e:
        raise SystemExit(e)


if __name__ == "__main__":
    main()
