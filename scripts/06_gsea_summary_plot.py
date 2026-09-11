import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

INPUT = Path(
    "results/gsea/custom_wnt_preranked_gsea.csv"
)

OUTPUT_DIR = Path("figures/gsea")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# Make labels presentation-friendly
labels = {
    "EPITHELIAL_JUNCTION_CORE": "Epithelial junction core",
    "TCF_Dependent": "TCF/LEF dependent",
    "BetaCatenin_Dependent": "β-catenin dependent",
    "BetaCatenin_Independent": "β-catenin independent",
    "TCF_Independent": "TCF/LEF independent",
    "WNT_PCP_CORE": "Wnt/PCP core",
}

df["label"] = df["Term"].map(labels).fillna(df["Term"])

df = df.sort_values("NES")

fig, ax = plt.subplots(figsize=(9, 5.5))
colors = ["#0072B2" if q < 0.05 else "#BDBDBD" for q in df["FDR q-val"]]

bars = ax.barh(
    df["label"],
    df["NES"],
    color=colors,
    edgecolor="black",
    linewidth=0.5,
)

# Mark FDR-significant result
for bar, q in zip(bars, df["FDR q-val"]):
    if q < 0.05:
        bar.set_hatch("//")

ax.axvline(
    0,
    linewidth=1,
    linestyle="--"
)

ax.set_xlabel("Normalized Enrichment Score (NES)")
ax.set_ylabel("")
ax.set_title(
    "Preranked GSEA of Epithelial, Canonical and Noncanonical Wnt-Related Gene Sets"
)

# Add q-values
for i, (_, row) in enumerate(df.iterrows()):

    x = row["NES"]

    offset = 0.04 if x >= 0 else -0.04
    ha = "left" if x >= 0 else "right"

    ax.text(
        x + offset,
        i,
        f"FDR={row['FDR q-val']:.3f}",
        va="center",
        ha=ha,
        fontsize=9
    )

ax.legend(
    handles=[
        Patch(facecolor="#0072B2", edgecolor="black", label="FDR q < 0.05"),
        Patch(facecolor="#BDBDBD", edgecolor="black", label="Not FDR-significant"),
    ],
    frameon=False,
    loc="lower right",
)
fig.text(0.5, 0.01, "Negative NES indicates enrichment toward genes decreased in mutant versus control.", ha="center", fontsize=9)

plt.tight_layout(rect=[0, 0.05, 1, 1])

png = OUTPUT_DIR / "wnt_gsea_NES_summary.png"
pdf = OUTPUT_DIR / "wnt_gsea_NES_summary.pdf"

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
