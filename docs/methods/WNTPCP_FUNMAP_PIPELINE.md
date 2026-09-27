# Auditable Wnt/PCP–FunMap overlay pipeline

This command-line workflow overlays **authoritative mouse DESeq2 results** on
**imported human FunMap functional-association scores** after an explicit
consensus one-to-one mouse–human orthology join. Optional STRING and curated
literature evidence are imported as separate layers and never relabeled as
FunMap evidence.

## Scientific safeguards

- The program does **not** infer coexpression or Mutual Rank from a DESeq2
  summary table. Those calculations require sample-level expression data.
- FunMap values are imported as predicted functional-association scores; they
  are not p-values and are not PP/PD/expression/activation annotations.
- Current-study tissue specificity belongs to the DESeq2 node attributes, not
  to the imported FunMap network.
- STRING confidence and literature claims remain separate from FunMap.
- Mouse and human symbols are never equated by capitalization. Missing
  orthologs are recorded; `--strict-orthology` makes them fatal.
- Gene-level BH adjustment is calculated over the complete valid DESeq2
  testing universe before display-oriented symbol deduplication.

## Required input schemas

### DESeq2

CSV or TSV containing at least:

```text
symbol,log2FoldChange,pvalue
```

`gene_id`, `baseMean`, and `padj` are optional. If `padj` is absent, the
pipeline calculates `padj_plot` using Benjamini–Hochberg correction across all
valid p-values before gene-symbol deduplication.

### Consensus one-to-one orthology

```text
human_symbol,mouse_symbol
```

The table must truly be one-to-one after case normalization. The program stops
if duplicated human or mouse symbols imply ambiguous mappings.

### FunMap

Preferred official format:

```text
gene1<TAB>gene2<TAB>score
```

The official FunMap repository describes this score as the predicted
probability of functional association. A documented derived table containing
`Gene1`, `Gene2`, and `Probability` is accepted but marked
`derived_funmap_probability_table` in provenance and edge outputs.

### Optional STRING

The loader accepts endpoint fields such as `preferredName_A` / `preferredName_B`
or `source` / `target`, with `score`, `combined_score`, or `string_score`.
Scores above 1 are divided by 1000. Use `--string-species human` when endpoints
are human; otherwise the default is mouse.

### Optional literature

See `literature_edges_template.csv`. `citation` should contain a DOI, PMID, or
stable URL. Literature rows are curated evidence and never converted into
FunMap probabilities.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Example

```bash
python wntpcp_funmap_pipeline.py \
  --deseq DESeq2_results.csv \
  --orthology consensus_mouse_human_one2one.csv \
  --funmap funmap.tsv \
  --string string_edges.tsv \
  --literature literature_edges.csv \
  --strict-orthology \
  --funmap-score-min 0.50 \
  --string-score-min 0.70 \
  --outdir results/wntpcp_funmap
```

For fish shell, the same command can be written with backslash continuations;
do not paste filenames or `.gitignore` patterns as standalone commands.

## Outputs

- `wntpcp_funmap_network.png` (600 dpi)
- `wntpcp_funmap_network.pdf` (vector)
- `wntpcp_funmap_network.svg` (editable vector)
- `wntpcp_funmap_network_edge_audit.csv`
- `wntpcp_funmap_network_node_audit.csv`
- `wntpcp_funmap_network_run_provenance.json`

The edge audit is the authoritative evidence record. Each row reports FunMap,
STRING, and literature fields separately, plus the current-study DESeq2
attributes of both endpoints. Figure labels use `F`, `S`, and `L` only as
compact layer codes.

## Interpretation

The figure is a **cross-species functional-association overlay**. It may support
hypothesis generation, but does not demonstrate direct physical binding,
regulator activity, edge direction, tissue-specific interaction, or causal
regulation in the current experiment.

The bundled `genes_default.txt` includes both the Wnt/PCP panel and the
structural/junction module (`Dsp`, `Pkp1`, `Dsg2`, `Cdh1`, `Jup`, `Krt5`,
`Krt14`, `Krt17`, `Bpifb3`, `Ano9`, `Aqp3`, `Il31ra`, `Gm6565`, and `Bmp6`).
Pass it with `--gene-file genes_default.txt`. The DESeq2 values are always read
from the authoritative table; the file contains symbols only.

## References

1. Zhang laboratory FunMap software and format documentation:
   <https://github.com/bzhanglab/funmap>
2. Zhang laboratory FunMap resource and downloads:
   <https://www.funmapnetwork.org/>
3. Dou Y. et al. *Mapping the functional network of human cancer through
   machine learning and pan-cancer proteogenomics.* Nature Cancer (2024).
   DOI: <https://doi.org/10.1038/s43018-024-00869-z>
4. STRING database: <https://string-db.org/>

## Important project-specific note

The supplied `human_mouse_consensus_one2one.csv` contains 162 mappings curated
for the Doumpas signature analysis and does not cover the default Wnt/PCP panel.
Use a broader audited one-to-one map for this network. Running with the narrow
file is still useful as a QC demonstration because the node audit will identify
every missing mapping; `--strict-orthology` will stop immediately.
