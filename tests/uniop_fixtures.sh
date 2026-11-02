#!/bin/sh
# Writes tests/data/uniop_example/ with the UniOP script's output on example/genome.gff.
# usage: UNIOP=/path/to/UniOP/src/UniOP [UNIOP_PYTHON=python] sh tests/uniop_fixtures.sh
# UniOP needs scikit-learn; the fork used is github.com/JulianTrouillon/UniOP (pair labels fixed for several sequences).
set -e
cd "$(dirname "$0")/.."
mkdir -p tests/data/uniop_example
"${UNIOP_PYTHON:-python}" "$UNIOP" -f example/genome.gff --gff_source ncbi_gff -t tests/data/uniop_example
