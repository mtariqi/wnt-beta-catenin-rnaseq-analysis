import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

INPUT = Path(
    "results/wnt_pcp_mechanistic/"
    "mechanistic_group_summary.csv"
)

OUTDIR = Path("figures/mechanistic")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

labels = {
    "EPITHELIAL": "Epithelial",
    "JUNCTION_ADHESION": "Junction / adhesion",
    "CYTOSKELETON_MIGRATION": "Cytoskeleton / migration",
    "STRESS_IFN": "Stress / interferon",
    "PCP_POLARITY": "PCP / polarity",
    "WNT_BETA_CATENIN": "Wnt / β-catenin",
}

df["Label"] = (
    df["Mechanistic_Group"]
    .map(labels)
    .fillna(df["Mechanistic_Group"])
)

# Order by fraction of significant pathways
df["Significant_Fraction"] = (
    df["significant_pathways"]
    / df["total_pathways"]
)

df = df.sort_values(
    "Significant_Fraction",
    ascending=True
)

fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.barh(
    df["Label"],
    df["Significant_Fraction"] * 100,
    color="#BDBDBD",
    edgecolor="black",
    linewidth=0.6,
)

ax.set_xlabel(
    "Pathways reaching FDR < 0.05 (%)"
)

ax.set_ylabel("")

ax.set_title(
    "FDR-Significant Pathways by Exploratory Mechanistic Group"
)

for bar, (_, row) in zip(
    bars,
    df.iterrows()
):

    width = bar.get_width()

    label = (
    f"{int(row['significant_pathways'])}/"
    f"{int(row['total_pathways'])} "
    f"({row['Significant_Fraction'] * 100:.1f}%)"
)
    ax.text(
        width + 0.3,
        bar.get_y() + bar.get_height()/2,
        label,
        va="center",
        fontsize=10
    )

ax.set_xlim(
    0,
    max(df["Significant_Fraction"] * 100) + 8,
)

fig.text(
    0.5,
    0.01,
    "No targeted pathway reached FDR q < 0.05; FDR was calculated within each source library.",
    ha="center",
    fontsize=9,
)

plt.tight_layout(rect=[0, 0.05, 1, 1])

png = (
    OUTDIR /
    "mechanistic_pathway_support.png"
)

pdf = (
    OUTDIR /
    "mechanistic_pathway_support.pdf"
)

plt.savefig(
    png,
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    pdf,
    bbox_inches="tight"
)

plt.close()

print(f"Saved: {png}")
print(f"Saved: {pdf}")
