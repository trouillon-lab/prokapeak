#!/usr/bin/env python3
"""Predict operons from gene spacing, following UniOP (Su, Zhang and Söding, bioRxiv 2024, 10.1101/2024.11.11.623000):
the intergenic distances of same-strand gene pairs are compared with those of convergent and divergent pairs, which
cannot be in one operon. Intended to reproduce the UniOP script's probabilities; checked in tests/."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from assign_targets import ext, read_genome

CUTOFF, N_NULL, SEED = 0.5, 10000, 0  # pair probability that joins two genes; null sample size and seed (as UniOP)
Q_RANGE = (0.01, 0.99)  # the prior (fraction of same-strand pairs in operons) is clipped to this range
MAX_H, MIN_OPPOSITE = 0.2, 20  # bandwidth above which probabilities tend to the prior; opposite-strand pairs for the null


def gene_pairs(g):
    """Consecutive genes along each chromosome (file order of chromosomes, then start):
    chrom, a, b, gap = start of b - end of a + 1 (as UniOP), class same / convergent / divergent."""
    g = g.assign(c=pd.factorize(g.chrom)[0]).sort_values(["c", "start"], kind="stable")
    a, b = g.iloc[:-1].reset_index(drop=True), g.iloc[1:].reset_index(drop=True)
    cls = np.where(a.strand == b.strand, "same", np.where(a.strand == "+", "convergent", "divergent"))
    p = pd.DataFrame({"chrom": a.chrom, "a": a.locus, "b": b.locus, "gap": b.start - a.end + 1, "class": cls})
    return p[a.chrom.to_numpy() == b.chrom.to_numpy()].reset_index(drop=True)


def prior(pairs):
    """q = (S - O) / S: same-strand pairs S minus opposite-strand pairs O, assuming same-strand pairs that are not in an
    operon are as frequent as opposite-strand ones. Returns the clipped value and the raw one."""
    s = (pairs["class"] == "same").sum()
    q = (2 * s - len(pairs)) / s if s else np.nan
    return min(max(q, Q_RANGE[0]), Q_RANGE[1]), q


def kde(x, sample, h):
    """Gaussian kernel density of sample, bandwidth h, at the points x (as scikit-learn's KernelDensity)."""
    out = np.empty(len(x))
    for i in range(0, len(x), 500):  # bounds the x by sample matrix
        d = (x[i:i + 500, None] - sample) / h
        out[i:i + 500] = np.exp(-0.5 * d * d).sum(1) / (len(sample) * h * np.sqrt(2 * np.pi))
    return out


def predict_operons(g, cutoff=CUTOFF, n_null=N_NULL, seed=SEED):
    """Gene table from read_genome -> operon table (operon id, comma-separated locus tags 5'->3'), pair table
    (gene_pairs columns plus prob, empty for opposite-strand pairs), list of warnings."""
    w, pairs = [], gene_pairs(g)
    cls, gap = pairs["class"].to_numpy(), pairs.gap.to_numpy()
    same, conv, div = (gap[cls == c] for c in ("same", "convergent", "divergent"))
    s, o = len(same), len(conv) + len(div)
    pairs["prob"] = np.nan
    if s == 0 or o == 0:
        w.append(f"predict: {s} same-strand and {o} opposite-strand pairs; no operon can be predicted")
        return pd.DataFrame(columns=["operon", "genes"]), pairs, w
    if o < MIN_OPPOSITE:
        w.append(f"predict: only {o} opposite-strand pairs to draw the null distances from")
    q, raw = prior(pairs)
    if q != raw:
        w.append(f"predict: prior {raw:.2f} clipped to {q}")
    np.random.seed(seed)  # legacy generator and call order as UniOP, so the null sample is the same
    null = [np.random.choice(x, n_null) for x in (conv, div) if len(x)]
    if len(null) < 2:
        w.append("predict: no convergent or no divergent pairs; the null distances come from the other class alone")
    null = np.mean(null, axis=0)  # mean of one convergent and one divergent distance
    srt = np.sort(same)
    rank = lambda x: (np.searchsorted(srt, x, side="right") + 0.5) / (s + 1)  # fraction of same-strand gaps <= x
    h = 2 * float(f"{100 / s:.2f}")  # UniOP's bandwidth
    if h == 0:
        h = 0.01
        w.append("predict: UniOP's bandwidth is 0 for this many pairs; 0.01 used")
    if h > MAX_H:
        w.append(f"predict: only {s} same-strand pairs; bandwidth {h:.2f} on the rank scale smooths the densities, "
                 f"probabilities tend to the prior {q:.2f}")
    r = rank(same)
    p = 1 - (1 - q) * kde(r, rank(null), h) / kde(r, r, h)  # 1 - P(pair is not in an operon | gap)
    p[p <= 0] = 0.01
    order = np.argsort(same, kind="stable")  # by gap: beyond the first minimum, p is held at that minimum (as UniOP)
    ps = p[order]
    ps[ps.argmin() + 1:] = ps.min()
    p[order] = ps
    pairs.loc[cls == "same", "prob"] = p

    join = pairs.prob.to_numpy() >= cutoff
    run = np.cumsum(~join | (pairs.chrom != pairs.chrom.shift()).to_numpy())  # consecutive joined pairs share a run
    strand = dict(zip(g.locus, g.strand))
    ops = []
    for _, u in pairs[join].groupby(run[join], sort=False):
        genes = u.a.tolist() + [u.b.iloc[-1]]
        ops.append(",".join(genes if strand[genes[0]] == "+" else genes[::-1]))
    return pd.DataFrame({"operon": [f"pred{i + 1}" for i in range(len(ops))], "genes": ops}), pairs, w


def main(argv=None):
    lst = lambda x: x.split(",")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff", required=True, help="GFF3 genome annotation")
    ap.add_argument("--gff-feature", type=lst, default="gene,CDS", help="GFF feature types to read")
    ap.add_argument("--gff-attribute", default="locus_tag", help="GFF attribute holding the gene id")
    ap.add_argument("--out", required=True, help="operon table (.tsv): operon id, comma-separated locus tags")
    ap.add_argument("--pairs", help="pair table (.tsv): chrom, gene a, gene b, gap, class, probability")
    ap.add_argument("--cutoff", type=float, default=CUTOFF, help="pair probability that joins two genes")
    a = ap.parse_args(argv)
    try:
        ext(a.gff, ".gff", ".gff3")
        ext(a.out, ".tsv")
        g, w = read_genome(a.gff, a.gff_feature, a.gff_attribute)
        ops, pairs, w2 = predict_operons(g, a.cutoff)
        for m in w + w2:
            print("warning:", m, file=sys.stderr)
        ops.to_csv(a.out, sep="\t", header=False, index=False)
        if a.pairs:
            pairs.to_csv(a.pairs, sep="\t", index=False, float_format="%.6f")
        q, _ = prior(pairs)
        print(f"{len(ops)} operons with {ops.genes.str.count(',').sum() + len(ops)} genes; "
              f"{(pairs['class'] == 'same').sum()} same-strand pairs, prior {q:.2f}, cutoff {a.cutoff}")
    except ValueError as e:
        raise SystemExit(e)


if __name__ == "__main__":
    main()
