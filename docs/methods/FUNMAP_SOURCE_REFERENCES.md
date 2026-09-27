# Source and evidence references

## FunMap

- Zhang laboratory software repository and output specification:
  <https://github.com/bzhanglab/funmap>
- Official FunMap network resource and downloadable edge scores:
  <https://www.funmapnetwork.org/>
- Dou Y, et al. Mapping the functional network of human cancer through machine
  learning and pan-cancer proteogenomics. *Nature Cancer*. 2024.
  <https://doi.org/10.1038/s43018-024-00869-z>

FunMap's official three-column edge output contains `gene1`, `gene2`, and a
predicted functional-association `score`. This pipeline retains that score as
`funmap_score`. It does not convert the score into a p-value or an experimental
interaction category.

## STRING

- STRING database: <https://string-db.org/>

STRING is an optional, independently labeled association layer. Its combined
confidence score is normalized to `[0,1]` when the imported values use the
usual `[0,1000]` scale. It is not merged numerically with the FunMap score.

## Orthology

The workflow requires a project-supplied, auditable consensus one-to-one table
with `human_symbol` and `mouse_symbol`. Mapping provenance should identify the
underlying releases or queries (for example, Ensembl and MGI) outside this
script. The loader validates one-to-one uniqueness but does not certify the
biological provenance of the supplied mapping.

## Current-study inference

DESeq2 log2 fold changes and p-values are read from the authoritative exported
table. When a source `padj` column is absent, `padj_plot` is generated with the
Benjamini–Hochberg procedure over all valid p-values before symbol
deduplication. That adjustment concerns current-study genes only and must not
be confused with FunMap or STRING scores.
