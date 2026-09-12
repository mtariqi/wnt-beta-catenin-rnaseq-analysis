#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "$SCRIPT_DIR/data" ]]; then
  PROJECT_ROOT="$SCRIPT_DIR"
else
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

DEFAULT_INPUT="$PROJECT_ROOT/data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
if [[ ! -f "$DEFAULT_INPUT" && -f "$PROJECT_ROOT/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM(1).csv" ]]; then
  DEFAULT_INPUT="$PROJECT_ROOT/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM(1).csv"
fi
INPUT="${1:-$DEFAULT_INPUT}"
OUTDIR="${2:-$PROJECT_ROOT/results/tf_sensitivity}"

python "$SCRIPT_DIR/14_upstream_tf_threshold_sensitivity.py" \
  --input "$INPUT" \
  --outdir "$OUTDIR" \
  --organism mouse \
  --enrichment-fdr 0.05 \
  --resume

python - <<'PY' "$OUTDIR"
import json
import sys
from pathlib import Path
import pandas as pd

p = Path(sys.argv[1])
counts = pd.read_csv(p / "threshold_input_counts.csv")
summary = pd.read_csv(p / "scenario_enrichment_summary.csv")
stability = pd.read_csv(p / "regulator_threshold_stability.csv")
agreement = pd.read_csv(p / "cross_database_agreement_by_scenario.csv")
metadata = json.loads((p / "run_metadata.json").read_text())

assert set(counts["scenario"]) == {"p005", "p010", "top2000", "top5000"}
assert counts.groupby("scenario")["direction"].nunique().eq(2).all()
assert counts.loc[counts.scenario.eq("top2000"), "gene_count"].sum() == 2000
assert counts.loc[counts.scenario.eq("top5000"), "gene_count"].sum() == 5000
assert summary["minimum_adjusted_p"].between(0, 1).all()
assert stability["scenarios_FDR_significant"].between(0, 4).all()
assert agreement["significant_in_both_databases"].isin([True, False]).all()

print("\nValidation passed.")
print("DESeq2 BH-audit significant symbols:", metadata["DESeq2_BH_FDR_lt_0.05_symbols"])
print("\nInput counts:\n", counts.to_string(index=False))
print("\nEnrichment summary:\n", summary.to_string(index=False))
print("\nMost stable regulator/library/direction combinations:\n",
      stability.head(25).to_string(index=False))
PY

echo
echo "Professor-response outputs are in: $OUTDIR"
