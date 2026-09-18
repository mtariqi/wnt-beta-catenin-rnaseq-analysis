# Professor-requested TF enrichment sensitivity analysis

This workflow compares four gene-list definitions without an absolute
log2-fold-change cutoff:

1. nominal p < 0.05;
2. nominal p < 0.10;
3. the 2,000 genes with the smallest nominal p-values; and
4. the 5,000 genes with the smallest nominal p-values.

The supplied dataset produces these input counts before enrichment:

| Scenario | Downregulated | Upregulated | Total |
|---|---:|---:|---:|
| Nominal p < 0.05 | 994 | 340 | 1,334 |
| Nominal p < 0.10 | 1,608 | 689 | 2,297 |
| Top 2,000 | 1,428 | 572 | 2,000 |
| Top 5,000 | 3,250 | 1,750 | 5,000 |

## Install

Place both scripts in the repository's `scripts/` directory, activate the
project virtual environment, and install the required packages:

```bash
python -m pip install pandas numpy matplotlib gseapy
chmod +x scripts/run_professor_tf_sensitivity.sh
```

## Audit the input lists without internet access

```bash
python scripts/14_upstream_tf_threshold_sensitivity.py \
  --input data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv \
  --outdir results/tf_sensitivity_audit \
  --audit-only
```

## Run the complete Enrichr sensitivity analysis

```bash
bash scripts/run_professor_tf_sensitivity.sh
```

The wrapper accepts optional input and output paths:

```bash
bash scripts/run_professor_tf_sensitivity.sh \
  data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv \
  results/tf_sensitivity
```

`--resume` is enabled by the wrapper. If a network request is interrupted,
rerunning the command reuses completed `*_all_results.csv` files.

## Main outputs

- `threshold_input_counts.csv`: exact input counts for every scenario.
- `threshold_input_counts.png`: input-count comparison figure.
- `combined_all_enrichment_results.csv`: complete Enrichr records.
- `combined_FDR05_enrichment_results.csv`: term-level BH-FDR results.
- `scenario_enrichment_summary.csv`: record counts and minimum FDR by scenario,
  direction, and database.
- `scenario_significant_record_counts.png`: comparison figure.
- `regulator_threshold_stability.csv`: recurrence across all four scenarios.
- `cross_database_agreement_by_scenario.csv`: TRRUST/ChEA replication audit.
- `run_metadata.json`: software, thresholds, and inferential-scope record.

## Interpretation boundary

Nominal p-values and top-N rankings define the input genes. The Enrichr
`Adjusted P-value` is a term-level multiple-testing correction for a particular
database, direction, and scenario. It must not be presented as a DESeq2 gene
FDR. Target-set enrichment also does not establish transcription-factor
activity, direct binding, regulatory direction, or causality.
