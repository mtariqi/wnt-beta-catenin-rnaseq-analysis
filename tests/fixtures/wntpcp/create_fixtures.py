#!/usr/bin/env python3
"""Generate non-sensitive synthetic inputs for workflow smoke tests."""

from pathlib import Path

import pandas as pd


OUTDIR = Path(__file__).resolve().parent
OUTDIR.mkdir(parents=True, exist_ok=True)


pd.DataFrame(
    {
        "gene_id": ["m1", "m2", "m3"],
        "symbol": ["Wnt5a", "Ror2", "Dsp"],
        "baseMean": [100, 80, 120],
        "log2FoldChange": [1.2, 0.8, -2.1],
        "pvalue": [0.01, 0.04, 0.001],
    }
).to_csv(OUTDIR / "deseq_example.csv", index=False)


pd.DataFrame(
    {
        "human_symbol": ["WNT5A", "ROR2", "DSP"],
        "mouse_symbol": ["Wnt5a", "Ror2", "Dsp"],
    }
).to_csv(OUTDIR / "orthology_example.csv", index=False)


pd.DataFrame(
    {
        "gene1": ["WNT5A", "WNT5A"],
        "gene2": ["ROR2", "DSP"],
        "score": [0.93, 0.61],
    }
).to_csv(
    OUTDIR / "funmap_example.tsv",
    sep="\t",
    index=False,
)


pd.DataFrame(
    {
        "preferredName_A": ["Wnt5a"],
        "preferredName_B": ["Dsp"],
        "combined_score": [850],
    }
).to_csv(
    OUTDIR / "string_example.tsv",
    sep="\t",
    index=False,
)


pd.DataFrame(
    {
        "source": ["Ror2"],
        "target": ["Dsp"],
        "species": ["mouse"],
        "relationship": ["curated association"],
        "citation": ["PMID:123"],
    }
).to_csv(
    OUTDIR / "literature_example.csv",
    index=False,
)


(OUTDIR / "genes_example.txt").write_text(
    "Wnt5a\nRor2\nDsp\n",
    encoding="utf-8",
)


print(f"Wrote synthetic workflow fixtures to {OUTDIR}")
