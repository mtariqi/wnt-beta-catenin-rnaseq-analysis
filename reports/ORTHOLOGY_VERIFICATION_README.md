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
- Consensus one-to-one mappings accepted: 162 (70.1%)
- Excluded or unresolved: 69
- Signature retention: 66.2%–71.9%
- Three genes (`FERMT3`, `PLAT`, and `CALCOCO1`) remained unresolved because Ensembl REST timed out; MGI-only assignments were deliberately not accepted.

## Validation result after verified mapping

All four signatures had negative normalized enrichment scores, consistent with a broad tendency toward lower expression in the current study. The effects were modest and not statistically significant after Benjamini–Hochberg correction:

| Signature | Detected genes | NES | Empirical two-sided P | BH FDR |
|---|---:|---:|---:|---:|
| TCF/LEF independent | 19 | -1.154 | 0.232 | 0.393 |
| β-catenin dependent | 115 | -1.179 | 0.234 | 0.393 |
| β-catenin independent | 42 | -0.853 | 0.377 | 0.393 |
| TCF/LEF dependent | 96 | -0.830 | 0.393 | 0.393 |

The strongest threshold-based exploratory signals were enrichment of TCF/LEF-dependent and β-catenin-dependent genes among the top 5,000 downregulated genes (nominal P = 0.013 and 0.019), but both had BH FDR = 0.449. Therefore, the defensible conclusion is a reproducible directional tendency, not formal statistical confirmation.

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
