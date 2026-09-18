#!/usr/bin/env python3
"""Create a publication-style TF target-set enrichment spider map.

The interactive figure is built with Plotly Express. The radial coordinate is
−log10(BH-adjusted enrichment p value); the orange ring marks FDR = 0.05.
Node/marker direction refers to the DEG list used for enrichment, not inferred
TF activation or inhibition.
"""

from pathlib import Path
import math

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


OUTPUT_DIR = Path("tf_spider_map_outputs")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Verified values from verified_slide10_plotted_values.csv.
VALUES = [
    ("SUZ12", 9.1e-17, "Downregulated DEG targets"),
    ("MTF2", 1.5e-06, "Downregulated DEG targets"),
    ("EZH2", 1.4e-05, "Downregulated DEG targets"),
    ("SOX2", 3.3e-05, "Downregulated DEG targets"),
    ("JARID2", 1.5e-04, "Downregulated DEG targets"),
    ("RING1B", 2.4e-04, "Downregulated DEG targets"),
    ("LEF1", 3.4e-04, "Downregulated DEG targets"),
    ("BMI1", 4.6e-04, "Downregulated DEG targets"),
    ("SP7", 1.7e-03, "Upregulated DEG targets"),
    ("NUCKS1", 4.1e-03, "Upregulated DEG targets"),
    ("EP300", 8.1e-03, "Upregulated DEG targets"),
    ("IRF8", 1.3e-02, "Upregulated DEG targets"),
    ("PPARG", 3.0e-02, "Upregulated DEG targets"),
]

COLORS = {
    "Downregulated DEG targets": "#0072B2",
    "Upregulated DEG targets": "#D55E00",
}
FDR_CUTOFF = 0.05
THRESHOLD = -math.log10(FDR_CUTOFF)

df = pd.DataFrame(VALUES, columns=["TF", "FDR", "Input direction"])
df["Evidence"] = -df["FDR"].map(math.log10)
df["FDR label"] = df["FDR"].map(lambda x: f"{x:.2e}")

# A single connected polygon preserves the ranking across all 13 TFs.
# Direction is then overlaid as color-coded markers.
fig = px.line_polar(
    df,
    r="Evidence",
    theta="TF",
    line_close=True,
    markers=False,
    template="simple_white",
    title="Transcription-factor target-set enrichment evidence",
)
fig.update_traces(
    line=dict(color="#4D4D4D", width=2),
    fill="toself",
    fillcolor="rgba(77,77,77,0.10)",
    name="Enrichment evidence",
    hoverinfo="skip",
    showlegend=True,
)

for direction, color in COLORS.items():
    sub = df[df["Input direction"] == direction]
    fig.add_trace(
        go.Scatterpolar(
            r=sub["Evidence"],
            theta=sub["TF"],
            mode="markers",
            name=direction,
            marker=dict(size=10, color=color, line=dict(color="white", width=1)),
            customdata=sub[["FDR", "Input direction"]].to_numpy(),
            hovertemplate=(
                "<b>%{theta}</b><br>BH-adjusted FDR=%{customdata[0]:.2e}"
                "<br>%{customdata[1]}<extra></extra>"
            ),
        )
    )

theta_closed = df["TF"].tolist() + [df["TF"].iloc[0]]
fig.add_trace(
    go.Scatterpolar(
        r=[THRESHOLD] * len(theta_closed),
        theta=theta_closed,
        mode="lines",
        name="FDR = 0.05",
        line=dict(color="#E69F00", width=2, dash="dash"),
        hoverinfo="skip",
    )
)

rmax = math.ceil(df["Evidence"].max() + 0.8)
fig.update_layout(
    font=dict(family="Arial", size=14, color="#1A1A1A"),
    title=dict(
        text=(
            "Transcription-factor target-set enrichment evidence"
            "<br><sup>Radial scale: −log₁₀(BH-adjusted enrichment p value); "
            "orange ring: FDR = 0.05</sup>"
        ),
        x=0.5,
        xanchor="center",
        font=dict(size=20),
    ),
    polar=dict(
        bgcolor="white",
        radialaxis=dict(
            range=[0, rmax],
            tickmode="array",
            tickvals=list(range(0, rmax + 1, 2)),
            gridcolor="#D9D9D9",
            linecolor="#999999",
            tickfont=dict(size=11),
            title=dict(text="−log₁₀(FDR)", font=dict(size=13)),
        ),
        angularaxis=dict(
            direction="clockwise",
            rotation=90,
            gridcolor="#E6E6E6",
            linecolor="#999999",
            tickfont=dict(size=13),
        ),
    ),
    legend=dict(
        orientation="h",
        yanchor="top",
        y=-0.10,
        xanchor="center",
        x=0.5,
        font=dict(size=12),
    ),
    margin=dict(l=90, r=90, t=110, b=145),
    width=1100,
    height=900,
    paper_bgcolor="white",
)
fig.add_annotation(
    text=(
        "Marker color denotes the DEG input direction used for enrichment; it does not "
        "establish TF activation, inhibition, direct binding, or causality."
    ),
    xref="paper",
    yref="paper",
    x=0.5,
    y=-0.19,
    xanchor="center",
    showarrow=False,
    font=dict(family="Arial", size=11, color="#595959"),
)

html_path = OUTPUT_DIR / "tf_enrichment_spider_map_interactive.html"
png_path = OUTPUT_DIR / "tf_enrichment_spider_map.png"
svg_path = OUTPUT_DIR / "tf_enrichment_spider_map.svg"
pdf_path = OUTPUT_DIR / "tf_enrichment_spider_map.pdf"
csv_path = OUTPUT_DIR / "tf_enrichment_spider_map_values.csv"

fig.write_html(html_path, include_plotlyjs="cdn")
df.to_csv(csv_path, index=False)

try:
    fig.write_image(png_path, width=1100, height=900, scale=3)
    fig.write_image(svg_path, width=1100, height=900)
    fig.write_image(pdf_path, width=1100, height=900)
except Exception as exc:
    # A headless environment may not provide the Chrome binary required by
    # Kaleido. Produce matching static publication files from the same table.
    import matplotlib.pyplot as plt
    import numpy as np

    labels = df["TF"].tolist()
    values = df["Evidence"].tolist()
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    closed_angles = angles + angles[:1]
    closed_values = values + values[:1]

    static_fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"polar": True})
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.plot(closed_angles, closed_values, color="#4D4D4D", linewidth=2)
    ax.fill(closed_angles, closed_values, color="#4D4D4D", alpha=0.10)
    ax.plot(
        closed_angles,
        [THRESHOLD] * len(closed_angles),
        color="#E69F00",
        linestyle="--",
        linewidth=2,
        label="FDR = 0.05",
    )
    for direction, color in COLORS.items():
        idx = [i for i, d in enumerate(df["Input direction"]) if d == direction]
        ax.scatter(
            [angles[i] for i in idx],
            [values[i] for i in idx],
            s=70,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=5,
            label=direction,
        )
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=11, fontfamily="DejaVu Sans")
    ax.set_ylim(0, rmax)
    ax.set_yticks(list(range(0, rmax + 1, 2)))
    ax.set_yticklabels([str(x) for x in range(0, rmax + 1, 2)], fontsize=9)
    ax.set_rlabel_position(78)
    ax.grid(color="#D9D9D9", linewidth=0.8)
    ax.spines["polar"].set_color("#999999")
    static_fig.suptitle(
        "Transcription-factor target-set enrichment evidence",
        fontsize=17,
        fontweight="bold",
        y=0.98,
    )
    ax.set_title(
        "Radial scale: −log₁₀(BH-adjusted enrichment p value); orange ring: FDR = 0.05",
        fontsize=11,
        pad=28,
    )
    handles, legend_labels = ax.get_legend_handles_labels()
    static_fig.legend(
        handles,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.075),
        ncol=3,
        frameon=False,
        fontsize=10,
    )
    static_fig.text(
        0.5,
        0.025,
        "Marker color denotes DEG input direction; enrichment does not establish TF activity, binding, or causality.",
        ha="center",
        fontsize=9,
        color="#595959",
    )
    static_fig.subplots_adjust(top=0.84, bottom=0.18, left=0.10, right=0.90)
    static_fig.savefig(png_path, dpi=300, facecolor="white")
    static_fig.savefig(svg_path, facecolor="white")
    static_fig.savefig(pdf_path, facecolor="white")
    plt.close(static_fig)
    print("Kaleido/Chrome unavailable; matching static files were rendered from the same values.")

print(f"Saved: {html_path}")
print(f"Saved: {csv_path}")
for path in (png_path, svg_path, pdf_path):
    if path.exists():
        print(f"Saved: {path}")
