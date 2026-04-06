#!/usr/bin/env bash
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
BIN_DIR="$ROOT_DIR/build/src/progs"
TREE_FILE="$ROOT_DIR/../gex_lineage_tests/input/quinn.90.var.nex"
OUT_DIR="$ROOT_DIR/../gex_lineage_tests/output"

SIM_LOW="$OUT_DIR/sim_test_low"
FIT_LOW="$OUT_DIR/sim_test_low_fit"
EVAL_LOW="$OUT_DIR/sim_test_low_eval"
SIM_MED="$OUT_DIR/sim_test_med"
FIT_MED="$OUT_DIR/sim_test_med_fit"
EVAL_MED="$OUT_DIR/sim_test_med_eval"
SIM_REPRO_A="$OUT_DIR/sim_test_repro_a"
SIM_REPRO_B="$OUT_DIR/sim_test_repro_b"
FILTER_SMOKE="$OUT_DIR/sim_test_filter_smoke"

get_metric() {
  awk -F '\t' -v key="$1" '$1 == key { print $2 }' "$2"
}

python_check() {
  python3 - "$@"
}

# Test 1: low-noise latent simulation.
# This is the "easy" recovery case: simulate a small latent model with strong
# lineage signal and low observation noise, run gexLineage end to end, then
# evaluate whether the fitted model recovers the planted structure well.
"$BIN_DIR/gexSimulateLatent" \
  --trees "$TREE_FILE" \
  --outprefix "$SIM_LOW" \
  --k 3 \
  --n-genes 60 \
  --sigma-obs 0.02 \
  --sigma-latent 2.0,1.0,0.5 \
  --seed 7

"$BIN_DIR/gexLineage" \
  --trees "$TREE_FILE" \
  --expr "$SIM_LOW.expr.tsv" \
  --outprefix "$FIT_LOW" \
  --n-filter-trees 3 \
  --n-model-trees 3 \
  --n-perms 100 \
  --pca-var-threshold 0.95 \
  --seed 7

"$BIN_DIR/gexEvaluateSimulation" \
  --truth-prefix "$SIM_LOW" \
  --fit-prefix "$FIT_LOW" \
  --outprefix "$EVAL_LOW"

# Test 2: moderate-noise latent simulation.
# This uses the same underlying latent setup but increases observation noise so
# recovery should become harder. We expect metrics to remain reasonable, but to
# degrade relative to the low-noise case rather than improving or staying flat.
"$BIN_DIR/gexSimulateLatent" \
  --trees "$TREE_FILE" \
  --outprefix "$SIM_MED" \
  --k 3 \
  --n-genes 60 \
  --sigma-obs 0.4 \
  --sigma-latent 2.0,1.0,0.5 \
  --seed 9

"$BIN_DIR/gexLineage" \
  --trees "$TREE_FILE" \
  --expr "$SIM_MED.expr.tsv" \
  --outprefix "$FIT_MED" \
  --n-filter-trees 3 \
  --n-model-trees 3 \
  --n-perms 100 \
  --pca-var-threshold 0.95 \
  --seed 9

"$BIN_DIR/gexEvaluateSimulation" \
  --truth-prefix "$SIM_MED" \
  --fit-prefix "$FIT_MED" \
  --outprefix "$EVAL_MED"


# Metric sanity checks.
# Enforce that the low-noise case has strong recovery, and that the
# moderate-noise case is measurably worse on representative recovery metrics.
python_check "$LOW_SUBSPACE" "$LOW_GENE_CORR" "$MED_SUBSPACE" "$MED_GENE_CORR" <<'PY'
import sys

low_subspace = float(sys.argv[1])
low_gene = float(sys.argv[2])
med_subspace = float(sys.argv[3])
med_gene = float(sys.argv[4])

if low_subspace < 0.95:
    raise SystemExit(f"low-noise subspace similarity too low: {low_subspace}")
if low_gene < 0.98:
    raise SystemExit(f"low-noise gene covariance correlation too low: {low_gene}")
if med_subspace >= low_subspace:
    raise SystemExit(
        f"moderate-noise subspace similarity did not degrade: low={low_subspace}, med={med_subspace}"
    )
if med_gene >= low_gene:
    raise SystemExit(
        f"moderate-noise gene covariance correlation did not degrade: low={low_gene}, med={med_gene}"
    )
PY

printf 'Simulation validation checks passed.\n'
