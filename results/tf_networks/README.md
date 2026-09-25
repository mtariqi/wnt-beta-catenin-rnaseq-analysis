# TF–Wnt–osteogenesis network analysis

All study-gene calculations were performed locally. TF edges denote membership in public TRRUST/ChEA target sets that pass FDR < 0.05 and overlap with nominal DEGs; they do not establish direct regulation in this tissue. Pathway membership came from GO Biological Process 2023 and WikiPathways 2024 Mouse.

| metric                                    |   value |
|:------------------------------------------|--------:|
| All tested unique genes                   |   21715 |
| Nominal p < 0.05                          |    1334 |
| Upregulated nominal DEGs                  |     340 |
| Downregulated nominal DEGs                |     994 |
| Wnt annotation genes in universe          |     412 |
| Osteogenesis annotation genes in universe |     319 |
| Wnt TF–DEG edges plotted                  |      48 |
| Osteogenesis TF–DEG edges plotted         |      23 |
| TFs shared between maps                   |       2 |

## Methodological basis and interpretation

This is a custom implementation for the present dataset, not a reproduction of a single published network-analysis pipeline. TF-target over-representation is evaluated with one-sided Fisher exact tests and Benjamini-Hochberg FDR correction. TF target sets are obtained from TRRUST and ChEA; pathway annotations are obtained from GO Biological Process and WikiPathways. Network construction and graph statistics use NetworkX. Node degree is encoded by node size. Edge width is intentionally constant and does not encode binding strength, confidence, causality, or interaction strength. Polycomb coloring is a prespecified visual annotation only.

## Suggested references

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate. Journal of the Royal Statistical Society Series B 57:289-300.
- Lachmann A, et al. 2010. ChEA: transcription factor regulation inferred from integrating genome-wide ChIP-X experiments. Bioinformatics 26:2438-2444.
- Han H, et al. 2018. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. Nucleic Acids Research 46:D380-D386.
- The Gene Ontology Consortium. Cite the current GO resource paper used by the manuscript.
- Agrawal A, et al. 2024. WikiPathways 2024: next generation pathway database. Nucleic Acids Research 52:D679-D689.
- Hagberg AA, Schult DA, Swart PJ. 2008. Exploring network structure, dynamics, and function using NetworkX. Proceedings of the 7th Python in Science Conference.
