#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from statsmodels.stats.multitest import multipletests

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DE_FILE = PROJECT_ROOT / "data" / "processed" / "DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
CRANIO_FILE = PROJECT_ROOT / "results" / "craniofacial_skeletal" / "04_craniofacial_skeletal_DESeq2_intersection.csv"
WNT_FILE = PROJECT_ROOT / "results" / "craniofacial_skeletal" / "07_canonical_Wnt_reference_genes.csv"
OUTDIR = PROJECT_ROOT / "figures" / "nature_volcano_craniofacial"
OUTDIR.mkdir(parents=True, exist_ok=True)

LOG2FC_CUTOFF = 1.0
Y_MODE = "minus_log10"   # change to "raw_p" if Professor wants raw P on y-axis
MAX_AUTO_LABELS = 10

FORCE_DEV_LABELS = {
    "Bmp6","Lama5","Tbx1","Pitx1","Col10a1","Sox10","Wnt6","Sox5"
}
WNT_LABELS = {"Axin2","Tcf7","Apcdd1","Dkk1","Nkd1","Sp5"}

COLOR_MAP = {
    "Not significant": "#C7C7C7",
    "0.05 ≤ p < 0.10": "#E69F00",
    "Nominal p < 0.05 only": "#56B4E9",
    "BH-FDR < 0.05": "#0072B2",
}
CATEGORY_ORDER = [
    "Not significant",
    "0.05 ≤ p < 0.10",
    "Nominal p < 0.05 only",
    "BH-FDR < 0.05",
]

df = pd.read_csv(DE_FILE)
df.columns = df.columns.str.strip()
df = df.dropna(subset=["symbol","log2FoldChange","pvalue"]).copy()
df["symbol"] = df["symbol"].astype(str).str.strip()
df["BH_FDR_AUDIT"] = multipletests(df["pvalue"].values, method="fdr_bh")[1]
df = (
    df.sort_values(["pvalue", "gene_id"], ascending=[True, True])
    .drop_duplicates("symbol", keep="first")
    .copy()
)

df["Significance"] = "Not significant"
df.loc[(df["pvalue"] >= 0.05) & (df["pvalue"] < 0.10), "Significance"] = "0.05 ≤ p < 0.10"
df.loc[(df["pvalue"] < 0.05) & (df["BH_FDR_AUDIT"] >= 0.05), "Significance"] = "Nominal p < 0.05 only"
df.loc[df["BH_FDR_AUDIT"] < 0.05, "Significance"] = "BH-FDR < 0.05"

df["pvalue_plot"] = df["pvalue"].clip(lower=1e-300)
if Y_MODE == "minus_log10":
    df["y_value"] = -np.log10(df["pvalue_plot"])
    y_title = "−log₁₀(P)"
else:
    df["y_value"] = df["pvalue_plot"]
    y_title = "P value"

cranio = pd.read_csv(CRANIO_FILE)
tier1 = cranio[cranio["BH_FDR_AUDIT"] < 0.05].copy()
tier2 = cranio[(cranio["pvalue"] < 0.05) & (cranio["BH_FDR_AUDIT"] >= 0.05)].copy()
tier2["label_score"] = (
    -np.log10(tier2["pvalue"].clip(lower=1e-300))
    + 0.35 * tier2["log2FoldChange"].abs()
    + 0.08 * tier2["GO_term_count"]
)
tier2 = tier2.sort_values("label_score", ascending=False)

developmental_labels = set(tier1["symbol"]) | set(tier2["symbol"].head(MAX_AUTO_LABELS)) | FORCE_DEV_LABELS
developmental_labels &= set(df["symbol"])
wnt_labels = WNT_LABELS & set(df.loc[df["pvalue"] < 0.05, "symbol"])

df["Label_group"] = "Other"
df.loc[df["symbol"].isin(developmental_labels), "Label_group"] = "Craniofacial / skeletal"
df.loc[df["symbol"].isin(wnt_labels), "Label_group"] = "Canonical Wnt reference"

fig = px.scatter(
    df,
    x="log2FoldChange",
    y="y_value",
    render_mode="svg",
    color="Significance",
    category_orders={"Significance": CATEGORY_ORDER},
    color_discrete_map=COLOR_MAP,
    hover_name="symbol",
    hover_data={
        "baseMean": ":.1f",
        "log2FoldChange": ":.3f",
        "pvalue": ":.3e",
        "BH_FDR_AUDIT": ":.3e",
        "Label_group": True,
        "y_value": False,
    },
    labels={
        "log2FoldChange": "log₂ fold change (mutant vs control)",
        "y_value": y_title,
        "Significance": "",
    },
)

for trace in fig.data:
    if trace.name == "Not significant":
        trace.update(marker=dict(size=4, opacity=0.25, line=dict(width=0)))
    elif trace.name == "0.05 ≤ p < 0.10":
        trace.update(marker=dict(size=4.3, opacity=0.50, line=dict(width=0)))
    else:
        trace.update(marker=dict(size=5, opacity=0.72, line=dict(width=0)))

# No fold-change cutoff: show only the zero-effect reference line.
fig.add_vline(x=0, line_width=0.8, line_color="black")

if Y_MODE == "minus_log10":
    fig.add_hline(y=-np.log10(0.05), line_width=0.8, line_dash="dash", line_color="black")
    fig.add_hline(y=-np.log10(0.10), line_width=0.8, line_dash="dot", line_color="black")
else:
    fig.add_hline(y=0.05, line_width=0.8, line_dash="dash", line_color="black")
    fig.add_hline(y=0.10, line_width=0.8, line_dash="dot", line_color="black")

dev = df[df["symbol"].isin(developmental_labels)].copy()
fig.add_trace(go.Scatter(
    x=dev["log2FoldChange"], y=dev["y_value"], mode="markers",
    name="Craniofacial / skeletal",
    marker=dict(size=9, color="white", line=dict(width=1.5, color="black"), symbol="circle"),
    text=dev["symbol"],
    customdata=np.stack([dev["pvalue"], dev["BH_FDR_AUDIT"]], axis=-1),
    hovertemplate="<b>%{text}</b><br>log2FC=%{x:.3f}<br>P=%{customdata[0]:.3e}<br>BH-FDR audit=%{customdata[1]:.3e}<extra></extra>"
))

wnt_plot = df[df["symbol"].isin(wnt_labels)].copy()
fig.add_trace(go.Scatter(
    x=wnt_plot["log2FoldChange"], y=wnt_plot["y_value"], mode="markers",
    name="Canonical Wnt reference",
    showlegend=not wnt_plot.empty,
    marker=dict(size=8, color="white", line=dict(width=1.3, color="black"), symbol="diamond"),
    text=wnt_plot["symbol"],
    customdata=np.stack([wnt_plot["pvalue"], wnt_plot["BH_FDR_AUDIT"]], axis=-1),
    hovertemplate="<b>%{text}</b><br>log2FC=%{x:.3f}<br>P=%{customdata[0]:.3e}<br>BH-FDR audit=%{customdata[1]:.3e}<extra></extra>"
))

offsets = {
    "Bmp6": (8,12), "Lama5": (-8,10), "Tbx1": (-8,12), "Pitx1": (-8,12),
    "Col10a1": (-8,-14), "Sox10": (-8,10), "Wnt6": (-8,-14), "Sox5": (8,12),
    "Axin2": (8,12), "Tcf7": (8,-14), "Apcdd1": (-8,12), "Dkk1": (8,12),
    "Nkd1": (-8,-14), "Sp5": (8,12),
}

for _, row in dev.iterrows():
    gene = row["symbol"]
    xshift, yshift = offsets.get(gene, (7 if row["log2FoldChange"] >= 0 else -7, 9))
    fig.add_annotation(
        x=row["log2FoldChange"], y=row["y_value"],
        text=(f"<b><i>{gene}</i></b>" if row["BH_FDR_AUDIT"] < 0.05 else f"<i>{gene}</i>"),
        showarrow=False, xshift=xshift, yshift=yshift,
        xanchor="left" if xshift >= 0 else "right",
        font=dict(family="Arial", size=10, color="black")
    )

for _, row in wnt_plot.iterrows():
    gene = row["symbol"]
    xshift, yshift = offsets.get(gene, (7,8))
    fig.add_annotation(
        x=row["log2FoldChange"], y=row["y_value"], text=f"<i>{gene}</i>",
        showarrow=False, xshift=xshift, yshift=yshift,
        xanchor="left" if xshift >= 0 else "right",
        font=dict(family="Arial", size=9, color="#555555")
    )

fig.update_layout(
    template="simple_white",
    width=1100, height=760,
    font=dict(family="Arial", size=11, color="black"),
    title=None, paper_bgcolor="white", plot_bgcolor="white",
    margin=dict(l=90, r=35, t=25, b=80),
    legend=dict(
        title="", x=0.99, y=0.99, xanchor="right", yanchor="top",
        bgcolor="rgba(255,255,255,0.90)", borderwidth=0,
        font=dict(family="Arial", size=10)
    ),
)
fig.update_xaxes(
    title="log₂ fold change (mutant vs control)",
    showline=True, linewidth=1, linecolor="black",
    ticks="outside", showgrid=False, zeroline=False
)
fig.update_yaxes(
    title=y_title,
    showline=True, linewidth=1, linecolor="black",
    ticks="outside", showgrid=False, zeroline=False
)

if Y_MODE == "raw_p":
    fig.update_yaxes(
        type="log", autorange="reversed",
        tickvals=[1,0.1,0.05,0.01,0.001,0.0001,0.00001,0.000001],
        ticktext=["1","0.1","0.05","0.01","0.001","0.0001","0.00001","0.000001"]
    )

html_file = OUTDIR / "Figure1_craniofacial_skeletal_volcano.html"
svg_file = OUTDIR / "Figure1_craniofacial_skeletal_volcano.svg"
pdf_file = OUTDIR / "Figure1_craniofacial_skeletal_volcano.pdf"
png_file = OUTDIR / "Figure1_craniofacial_skeletal_volcano.png"

fig.write_html(html_file, include_plotlyjs="cdn")

label_table = pd.concat([
    dev.assign(Figure_group="Craniofacial / skeletal"),
    wnt_plot.assign(Figure_group="Canonical Wnt reference")
], ignore_index=True)

label_table[
    ["symbol","Figure_group","baseMean","log2FoldChange","pvalue","BH_FDR_AUDIT","Significance"]
].sort_values(["Figure_group","pvalue"]).to_csv(
    OUTDIR / "Figure1_labeled_gene_table.csv", index=False
)

try:
    fig.write_image(svg_file, width=1100, height=760)
    fig.write_image(pdf_file, width=1100, height=760)
    fig.write_image(png_file, width=1100, height=760, scale=3)
    print(f"Saved: {svg_file}")
    print(f"Saved: {pdf_file}")
    print(f"Saved: {png_file}")
except Exception as exc:
    print("\nStatic export failed. Install/update Kaleido:\npython -m pip install -U kaleido\n")
    print(exc)

print(f"Saved: {html_file}")

print("\nCraniofacial/skeletal genes highlighted:")
print(dev[["symbol","log2FoldChange","pvalue","BH_FDR_AUDIT","Significance"]]
      .sort_values(["BH_FDR_AUDIT","pvalue"]).to_string(index=False))

print("\nCanonical Wnt reference genes:")
print(wnt_plot[["symbol","log2FoldChange","pvalue","BH_FDR_AUDIT","Significance"]]
      .sort_values("pvalue").to_string(index=False))
