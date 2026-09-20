# Wnt paper-signature ranked validation

## Scope

This analysis follows Professor Wei Hsu's instruction to continue using the
published gene signatures as reference sets. It does **not** reanalyse the
paper's raw RNA-seq data.

The reference workbook supplies four Doumpas et al. (2019) gene sets and the
published wild-type (WT) log-fold changes and P values:

- beta-catenin independent
- beta-catenin dependent
- TCF/LEF independent
- TCF/LEF dependent

## Methods

1. The current 4-control versus 4-mutant DESeq2 summary table was filtered to
   21,715 genes with a valid symbol, base mean, log2 fold change, and nominal
   P value. Duplicate symbols retained the row with the smallest P value.
2. Genes were ranked by `sign(log2FoldChange) * -log10(P)`. Infinite scores
   were prevented by clipping zero P values to the smallest positive floating-
   point value and capping `-log10(P)` at 300.
3. Each paper signature was tested using 10,000 Monte Carlo random sets that
   preserved its distribution across current-study baseMean deciles. The
   normalized enrichment score (NES) is the observed mean rank score minus the
   null mean, divided by the null standard deviation.
4. Directional overlaps were tested for six selection scenarios: P<0.01,
   P<0.05, P<0.10, top 1,000, top 2,000, and top 5,000. Each was split into
   upregulated and downregulated genes. Fisher's exact test supplied the odds
   ratio and one-sided enrichment P value.
5. Direction concordance compared the sign of the paper's WT logFC with the
   sign of the current log2 fold change. A two-sided binomial test evaluated
   concordance against 50%; Spearman correlation evaluated continuous fold-
   change agreement.
6. Benjamini-Hochberg adjustment was applied within each output family.

## Results

- All four expression-matched ranked NES values were negative, ranging from
  -1.21 to -0.77, but none was significant after BH correction (minimum FDR
  0.335).
- The strongest threshold-specific nominal signals occurred in the top-5,000
  downregulated pool: TCF/LEF-dependent (OR 1.80; P=0.011) and beta-catenin-
  dependent (OR 1.65; P=0.018). Neither survived correction across the 48
  threshold/direction/signature tests (both FDR 0.429).
- Direction concordance was 59.2%-60.0% for all four signatures. The nominal
  binomial P value for the beta-catenin-dependent set was 0.035, but its BH FDR
  was 0.114. No Spearman correlation or concordance test was FDR significant.

## Interpretation

The results show a consistent but statistically modest tendency toward
downregulation and approximately 60% cross-study direction agreement. They do
not establish significant genome-wide enrichment after multiple-testing
correction. The top-5,000 findings are hypothesis-generating sensitivity
signals, not confirmed pathway effects. This transparent result supports
continuing with signature-based biological interpretation and targeted review;
it does not yet justify raw-paper-data reanalysis unless Professor Hsu decides
the present information is insufficient.

## Cross-species limitation

Without a curated ortholog table, uppercase symbol matches are labeled
"symbol-compatible; orthology not formally verified." The script accepts an
optional `--ortholog-map` CSV/TSV containing `human_symbol` and `mouse_symbol`.
One-to-many mappings are excluded automatically because they require manual
biological adjudication.

## Reproduction

```bash
python scripts/21_wnt_signature_ranked_validation.py \
  --deseq2 data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv \
  --signatures data/reference/Wnt_beta_catenin_LEFTCF_signatures.xlsx \
  --outdir results/wnt_signature_ranked_validation \
  --n-sim 10000 \
  --seed 20260917
```

## Reference

Doumpas, N., Lampart, F., Robinson, M. D., Lentini, A., Nestor, C. E., Cantù,
C., & Basler, K. (2019). TCF/LEF dependent and independent transcriptional
regulation of Wnt/beta-catenin target genes. *The EMBO Journal, 38*, e98873.
https://doi.org/10.15252/embj.201798873
