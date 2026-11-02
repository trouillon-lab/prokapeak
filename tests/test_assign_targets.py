"""assign_targets on example/: the expected outputs, the predicted-operon files, and the input checks."""
import filecmp
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from assign_targets import assign, build_units, main, read_genome, resolve

ROOT = Path(__file__).resolve().parents[1]
E = ROOT / "example"
PEAKS = sorted(str(p) for p in (E / "peaks").glob("*.narrowPeak"))


def test_example_expected(tmp_path, capsys):
    main(["--gff", str(E / "genome.gff"), "--operons", str(E / "operons.tsv"), "--tss", str(E / "tss.bed"),
          "--peaks", *PEAKS, "--out", str(tmp_path)])
    assert capsys.readouterr().err == ""
    for f in (E / "expected").glob("*.bed"):
        assert filecmp.cmp(f, tmp_path / f.name, shallow=False), f.name


def test_predict_operons_files(tmp_path):
    main(["--gff", str(E / "genome.gff"), "--predict-operons", "--peaks", PEAKS[0], "--out", str(tmp_path)])
    for name in ("predicted_operons.tsv", "predicted_operon_pairs.tsv"):
        assert filecmp.cmp(E / "expected" / name, tmp_path / name, shallow=False), name


@pytest.mark.parametrize("args, msg", [
    (["--operons", "wide.bed"], "must end in .tsv"),
    (["--operons", "missing.tsv"], "file not found"),
    (["--operons", str(E / "operons.tsv"), "--predict-operons"], "cannot be combined"),
    (["--peaks", "wide.bed"], "not 1-bp summit"),
])
def test_input_errors(tmp_path, args, msg):
    (tmp_path / "wide.bed").write_text("chrom1\t100\t200\n")
    args = [str(tmp_path / a) if a.endswith(".bed") or a == "missing.tsv" else a for a in args]
    if "--peaks" not in args:
        args += ["--peaks", PEAKS[0]]
    with pytest.raises(SystemExit, match=msg):
        main(["--gff", str(E / "genome.gff"), "--out", str(tmp_path)] + args)


def test_one_row_per_peak_and_gene():
    """Operon A,B,C with a 100-bp first gene: a summit 1 bp upstream of A is also in B's internal window."""
    gff = "".join(f"c\tt\tgene\t{a}\t{b}\t.\t+\t.\tlocus_tag={x}\n" for x, a, b in
                  [("A", 1001, 1100), ("B", 1121, 2000), ("C", 2021, 3000)])
    g, _ = read_genome(io.StringIO(gff))
    units, _ = build_units(g, pd.DataFrame([("op", "A,B,C")]))
    pk = pd.DataFrame([["c", 999, 1000, "pk", 0, ".", 0.0, 0.0, 0.0, 0]])
    bed, _, _ = assign(units, g, pk, np.array([1000]), "both", "closest")
    assert bed.target.tolist() == ["A", "B", "C"] and (bed.distance_to_anchor == -1).all()


@pytest.mark.parametrize("rule", ["closest", "farthest"])
def test_overlap_ties_keep_both(rule):
    """A summit midway between two divergent anchors: the order of the units must not decide."""
    h = pd.DataFrame({"p": 0, "u": [0, 1], "d": [150, 150], "s": [1, -1], "uid": [0, 1]})
    assert resolve(h, rule, rule).u.tolist() == [0, 1]


def test_overlap_rule_checked():
    h = pd.DataFrame({"p": 0, "u": [0, 1], "d": [10, 300], "s": 1, "uid": [0, 1]})
    with pytest.raises(ValueError, match="nearest"):
        resolve(h, "nearest")


# One summit (p = 0) in the windows of four units: + strand at 10 and 300 bp, - strand at 50 and 250 bp from their
# anchors. The units kept for each --same-strand / --opp-strand pair.
HITS = pd.DataFrame({"p": 0, "u": [0, 1, 2, 3], "d": [10, 300, 50, 250], "s": [1, 1, -1, -1], "uid": [0, 1, 2, 3]})


@pytest.mark.parametrize("same, opp, kept", [
    ("closest", "closest", [0]), ("closest", "both", [0, 2]), ("closest", "farthest", [2]),
    ("both", "closest", [0, 1]), ("both", "both", [0, 1, 2, 3]), ("both", "farthest", [2, 3]),
    ("farthest", "closest", [3]), ("farthest", "both", [1, 3]), ("farthest", "farthest", [1]),
])
def test_overlap_rules(same, opp, kept):
    assert resolve(HITS, same, opp).u.tolist() == kept
