# Formal human–mouse orthology verification for the Wnt signature analysis

## Purpose

This workflow replaces capitalization-based human-to-mouse symbol conversion with an auditable, conservative cross-species mapping. It uses the published human Wnt/β-catenin signatures as reference sets and the current mouse RNA-seq DESeq2 results as the test dataset. It does not reanalyse the paper's raw RNA-seq data.

## Acceptance rule

A human gene is accepted into the primary map only when all three conditions hold:

1. Ensembl Compara reports exactly one mouse `ortholog_one2one` target.
2. The MGI homology class contains exactly one human gene and one mouse gene.
3. The Ensembl and MGI mouse symbols agree.

One-to-many, many-to-many, missing, discordant, and query-unresolved mappings are excluded from the primary analysis and retained in the audit tables.

## Final mapping outcome

- Unique human genes evaluated: 231
- Consensus one-to-one mappings accepted: **163** (70.6%)
- Excluded or unresolved: 68
- Signature retention: 66.2%–72.7%
- Two genes (`PLAT` and `CALCOCO1`) remained query-unresolved because Ensembl REST timed out; MGI-only assignments were conservatively excluded.
- `FERMT3`, which was unresolved in the earlier run, was recovered after a successful Ensembl retry and passed the Ensembl–MGI consensus one-to-one acceptance rule.

## Validation result after verified mapping

All four signatures had negative normalized enrichment scores, indicating a shared directional tendency toward lower expression in the current study. However, none was statistically significant after Benjamini–Hochberg correction (minimum empirical two-sided P = 0.233; minimum FDR = 0.491). Therefore, these findings provide exploratory directional evidence rather than statistically confirmed signature enrichment.

| Signature | Detected genes | NES | Empirical two-sided P | BH FDR |
|---|---:|---:|---:|---:|
| TCF/LEF independent | 19 | -1.149 | 0.233 | 0.491 |
| β-catenin dependent | 116 | -1.090 | 0.276 | 0.491 |
| β-catenin independent | 42 | -0.853 | 0.368 | 0.491 |
| TCF/LEF dependent | 97 | -0.660 | 0.502 | 0.502 |

## Reproduction

```bash
python scripts/22_verify_human_mouse_orthology.py \
  --signatures data/reference/Wnt_beta_catenin_LEFTCF_signatures.xlsx \
  --outdir results/wnt_orthology_verification \
  --workers 5 \
  --timeout 45 \
  --retries 3

python scripts/21_wnt_signature_ranked_validation.py \
  --deseq2 data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv \
  --signatures data/reference/Wnt_beta_catenin_LEFTCF_signatures.xlsx \
  --ortholog-map results/wnt_orthology_verification/human_mouse_consensus_one2one.csv \
  --outdir results/wnt_signature_ranked_validation_verified_orthologs \
  --n-sim 10000 \
  --seed 20260917
```

The workbook warning from `openpyxl` concerns an unsupported Excel extension and does not prevent the gene tables from being read.

## Key files

- `human_mouse_consensus_one2one.csv`: primary accepted mapping.
- `wnt_signature_orthology_audit.csv`: full per-signature audit trail.
- `ensembl_all_candidate_orthologs.csv`: all returned Ensembl candidates.
- `excluded_or_ambiguous_mappings.csv`: every nonaccepted mapping and reason.
- `orthology_retention_by_signature.csv`: retention counts and fractions.
- `orthology_provenance.json`: URLs, rule, timestamp, and MGI file checksum.
- `ranked_signature_monte_carlo.csv`: ranked Monte Carlo results.
- `threshold_sensitivity.csv`: threshold sensitivity analysis.
- `cross_study_direction_concordance.csv`: direction and rank concordance.

## Integrity statement

The code and narrative were developed with LLM assistance. The investigator remains responsible for validating inputs, inspecting mappings and exclusions, executing the analysis, interpreting the biology, and reporting the results. No LLM-generated claim should replace review of the source paper, database documentation, or primary output tables.
