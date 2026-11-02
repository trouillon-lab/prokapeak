"""predict_operons against the UniOP script (fixtures in data/uniop_example, made by uniop_fixtures.sh; a live run
when UNIOP names the script and UNIOP_PYTHON an interpreter with scikit-learn) and on a synthetic genome."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from assign_targets import read_genome
from predict_operons import main, predict_operons, prior

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "example" / "genome.gff"
TOL = 1e-6  # UniOP prints six decimals
UNIOP_FEATURES = ("gene", "pseudogene")  # what UniOP reads with --gff_source ncbi_gff


def read_uniop(d):
    """UniOP output directory -> pair table (a, b, p by locus tag) and the set of operons (sorted locus tag tuples)."""
    genes = pd.read_csv(d / "uniop.genes", sep="\t", dtype=str)
    tag = dict(zip(genes.gene.astype(int), genes.locus_tag))
    pred = pd.read_csv(d / "uniop.pred", sep="\t", skiprows=1, dtype=str).set_axis(["a", "b", "p"], axis=1)
    pred["p"] = pd.to_numeric(pred.p.str.strip(), errors="coerce")  # blank = opposite-strand pair
    pred[["a", "b"]] = pred[["a", "b"]].astype(int).apply(lambda c: c.map(tag))
    ops = pd.read_csv(d / "uniop.operon").idx_genes
    return pred, {tuple(sorted(tag[int(x)] for x in s.strip("[]").split(","))) for s in ops}


def check_against_uniop(gff, d):
    g, _ = read_genome(gff, UNIOP_FEATURES)
    ops, pairs, _ = predict_operons(g)
    pred, u_ops = read_uniop(d)
    m = pairs.merge(pred, on=["a", "b"], how="inner")
    assert len(m) == len(pairs), "every consecutive pair within a chromosome has a UniOP row"
    assert (m["class"].eq("same") == m.p.notna()).all(), "same-strand pairs are exactly the scored rows"
    same = m[m["class"] == "same"]
    assert len(same) > 0 and (same.prob - same.p).abs().max() <= TOL
    assert {tuple(sorted(x.split(","))) for x in ops.genes} == u_ops


def test_matches_uniop_fixtures():
    check_against_uniop(EXAMPLE, ROOT / "tests" / "data" / "uniop_example")


@pytest.mark.skipif("UNIOP" not in os.environ, reason="set UNIOP=/path/to/UniOP/src/UniOP and UNIOP_PYTHON")
def test_matches_uniop_live(tmp_path):
    gff = EXAMPLE
    cmd = [os.environ.get("UNIOP_PYTHON", sys.executable), os.environ["UNIOP"], "-f", str(gff),
           "--gff_source", "ncbi_gff", "-t", str(tmp_path)]
    subprocess.run(cmd, check=True, capture_output=True)
    check_against_uniop(gff, tmp_path)


def synthetic(n_operons=1500, seed=1):
    """One chromosome of operons of 1 to 5 genes; short gaps inside operons, long gaps between them.
    Returns the gene table and the set of planted same-strand pairs (a, b)."""
    rng = np.random.default_rng(seed)
    rows, planted, pos = [], set(), 1000
    for i in range(n_operons):
        strand, size = rng.choice(["+", "-"]), rng.integers(1, 6)
        for j in range(size):
            length = rng.integers(300, 1500)
            rows.append((f"G{len(rows):04d}", "chr", pos, pos + length, strand))
            if j:
                planted.add((rows[-2][0], rows[-1][0]))
            pos += length + (rng.integers(-10, 30) if j < size - 1 else rng.integers(100, 400))
    return pd.DataFrame(rows, columns=["locus", "chrom", "start", "end", "strand"]), planted


def test_synthetic_planted_operons():
    g, planted = synthetic()
    ops, pairs, w = predict_operons(g)
    assert not w
    same = pairs[pairs["class"] == "same"]
    called = {(a, b) for a, b, p in zip(same.a, same.b, same.prob) if p >= 0.5}
    hit = len(called & planted)
    assert hit / len(called) >= 0.95 and hit / len(planted) >= 0.95  # precision, recall
    _, q = prior(pairs)
    assert abs(q - len(planted) / len(same)) < 0.1
    assert all(len(x.split(",")) >= 2 for x in ops.genes)


def test_no_opposite_strand_pairs():
    g, _ = synthetic(20)
    g["strand"] = "+"
    ops, pairs, w = predict_operons(g)
    assert ops.empty and pairs.prob.isna().all() and any("no operon" in m for m in w)


def test_no_pair_across_chromosomes():
    g, _ = synthetic(40)
    g.loc[g.index >= len(g) // 2, "chrom"] = "plasmid"
    _, pairs, _ = predict_operons(g)
    assert len(pairs) == len(g) - 2 and (pairs.chrom.values[:-1] <= pairs.chrom.values[1:]).all()
    assert set(zip(pairs.a, pairs.b)).isdisjoint({(g.locus[len(g) // 2 - 1], g.locus[len(g) // 2])})


def test_three_gene_genome_runs():
    g = pd.DataFrame({"locus": list("abc"), "chrom": "chr", "start": [1, 500, 1000], "end": [400, 900, 1400],
                      "strand": ["+", "+", "-"]})
    ops, pairs, w = predict_operons(g)
    assert len(pairs) == 2 and len(w) >= 1


def test_command_line_new_files(tmp_path):
    """The command writes new output files that match example/expected/."""
    main(["--gff", str(EXAMPLE), "--out", str(tmp_path / "ops.tsv"), "--pairs", str(tmp_path / "pairs.tsv")])
    for out, name in (("ops.tsv", "predicted_operons.tsv"), ("pairs.tsv", "predicted_operon_pairs.tsv")):
        assert (tmp_path / out).read_text() == (ROOT / "example" / "expected" / name).read_text(), name
