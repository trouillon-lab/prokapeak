"""assign_targets on example/: the expected outputs, the predicted-operon files, and the input checks."""
import filecmp
from pathlib import Path

import pytest

from assign_targets import main

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
    (["--operons", "x.bed"], "must end in .tsv"),
    (["--operons", str(E / "operons.tsv"), "--predict-operons"], "cannot be combined"),
    (["--peaks", "wide.bed"], "not 1-bp summit"),
])
def test_input_errors(tmp_path, args, msg):
    (tmp_path / "wide.bed").write_text("chrom1\t100\t200\n")
    args = [str(tmp_path / a) if a.endswith("wide.bed") else a for a in args]
    if "--peaks" not in args:
        args += ["--peaks", PEAKS[0]]
    with pytest.raises(SystemExit, match=msg):
        main(["--gff", str(E / "genome.gff"), "--out", str(tmp_path)] + args)
