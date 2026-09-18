#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS="${1:-$PROJECT_ROOT/results/tf_sensitivity}"
REPORT="${2:-$PROJECT_ROOT/reports/manuscripts/Wnt_beta_catenin_RNAseq_Graduate_Report_REVISED.docx}"

python "$SCRIPT_DIR/15_tf_threshold_stability_heatmap.py" \
  --input "$RESULTS/combined_all_enrichment_results.csv" \
  --outdir "$RESULTS"

python "$SCRIPT_DIR/16_append_tf_sensitivity_to_report.py" \
  --report "$REPORT" \
  --results "$RESULTS" \
  --output-report "$PROJECT_ROOT/reports/manuscripts/Wnt_TF_sensitivity_report_updated.docx" \
  --email "$PROJECT_ROOT/reports/Professor_Hsu_TF_sensitivity_email.txt"

echo "Reporting workflow completed."
