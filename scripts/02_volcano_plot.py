import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

INPUT = Path("data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv")
OUTPUT = Path("figures/volcano/volcano_mutant_vs_control.png")

P_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
df.columns = df.columns.str.strip()

df = df.dropna(subset=["symbol", "log2FoldChange", "pvalue"]).copy()

df["status"] = "Not significant"

df.loc[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"] >= LOG2FC_CUTOFF),
    "status"
] = "Upregulated"

df.loc[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"] <= -LOG2FC_CUTOFF),
    "status"
] = "Downregulated"

# Avoid log10(0)
df["minus_log10_p"] = -np.log10(
    df["pvalue"].clip(lower=np.finfo(float).tiny)
)

fig, ax = plt.subplots(figsize=(9, 7))

groups = [
    ("Not significant", 8, 0.25),
    ("Upregulated", 18, 0.75),
    ("Downregulated", 18, 0.75),
]

for group, size, alpha in groups:
    subset = df[df["status"] == group]

    ax.scatter(
        subset["log2FoldChange"],
        subset["minus_log10_p"],
        s=size,
        alpha=alpha,
        label=f"{group} (n={len(subset)})"
    )

ax.axvline(
    LOG2FC_CUTOFF,
    linestyle="--",
    linewidth=1
)

ax.axvline(
    -LOG2FC_CUTOFF,
    linestyle="--",
    linewidth=1
)

ax.axhline(
    -np.log10(P_CUTOFF),
    linestyle="--",
    linewidth=1
)

# Label strongest DEGs
significant = df[df["status"] != "Not significant"].copy()

top_genes = significant.nsmallest(15, "pvalue")

for _, row in top_genes.iterrows():
    ax.annotate(
        row["symbol"],
        (
            row["log2FoldChange"],
            row["minus_log10_p"]
        ),
        xytext=(3, 3),
        textcoords="offset points",
        fontsize=8
    )

ax.set_xlabel("log2 Fold Change (Mutant vs Control)", fontsize=11)
ax.set_ylabel("-log10(p-value)", fontsize=11)

ax.set_title(
    "Differential Gene Expression\n4 Mutant vs 4 Control",
    fontsize=13
)

ax.legend(
    frameon=False,
    fontsize=9
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    OUTPUT,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Volcano plot saved to: {OUTPUT}")
print()
print(df["status"].value_counts())
