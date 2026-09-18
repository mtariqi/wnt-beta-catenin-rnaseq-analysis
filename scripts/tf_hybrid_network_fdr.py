#!/usr/bin/env python3
"""
tf_hybrid_network_fdr.py

Combines the two previous figures into one: the literature-based network
of TF-Wnt-osteogenesis relationships (tf_wnt_osteogenesis_network.py),
with node SIZE and COLOR INTENSITY now encoding each TF's real,
already-verified BH-adjusted FDR from the Enrichr/GSEApy audit (the same
13 values used in tf_enrichment_spider_chart.py).

What each visual channel means:
    - Node position / which hub(s) it connects to = literature-based
      relationship (see tf_wnt_osteogenesis_network.py for the citation
      basis of each edge -- unchanged from that script).
    - Edge style (solid/dashed/dotted/antagonist-colored) = strength and
      nature of that literature relationship, exactly as before.
    - Node SIZE and color SATURATION = statistical significance
      (-log10 of the real BH-adjusted FDR), scaled globally across all 13
      TFs with a sqrt transform so the single most extreme value (SUZ12)
      doesn't visually flatten the other 12 into indistinguishable dots.

This means a node can now tell you two things at a glance: which pathway
the literature associates it with (hub position, edge style), and how
strong the actual enrichment evidence for that TF was in this dataset
(size/color). Neither channel is estimated -- edges are literature-based
(see the other script's docstring) and sizes are computed directly from
the verified FDR values below.

Dependencies:
    pip install plotly kaleido

Usage:
    python tf_hybrid_network_fdr.py
"""

import math

import plotly.graph_objects as go

# ---------------------------------------------------------------
# Okabe-Ito colorblind-safe palette (Wong 2011, Nature Methods)
# ---------------------------------------------------------------
COL_WNT_HUB = "#0072B2"
COL_OSTEO_HUB = "#D55E00"
COL_NODE_BORDER = "#333333"
COL_EDGE_STRONG = "#333333"
COL_EDGE_INDIRECT = "#8C8C8C"
COL_EDGE_WEAK = "#B3B3B3"
COL_EDGE_ANTAGONIST = "#CC79A7"

# ---------------------------------------------------------------
# Real, verified FDR values (same 13 as the enrichment bar chart /
# spider chart -- see verified_slide10_plotted_values.csv)
# ---------------------------------------------------------------
FDR = {
    "SUZ12": 9.1e-17, "MTF2": 1.5e-06, "EZH2": 1.4e-05, "SOX2": 3.3e-05,
    "JARID2": 1.5e-04, "RING1B": 2.4e-04, "LEF1": 3.4e-04, "BMI1": 4.6e-04,
    "SP7": 1.7e-03, "NUCKS1": 4.1e-03, "EP300": 8.1e-03, "IRF8": 1.3e-02, "PPARG": 3.0e-02,
}
DOWN_GROUP = {"SUZ12", "MTF2", "EZH2", "SOX2", "JARID2", "RING1B", "LEF1", "BMI1"}
UP_GROUP = {"SP7", "NUCKS1", "EP300", "IRF8", "PPARG"}

_neglog = {k: -math.log10(v) for k, v in FDR.items()}
_sqrt = {k: math.sqrt(v) for k, v in _neglog.items()}
_lo, _hi = min(_sqrt.values()), max(_sqrt.values())


def node_size(tf):
    """Global sqrt-scaled size, range ~18-58px, so SUZ12 stands out
    without visually erasing the other 12 TFs' real size differences."""
    norm = (_sqrt[tf] - _lo) / (_hi - _lo)
    return 18 + norm * 40


def node_color(tf):
    """Color intensity within the TF's own group hue (blue for the
    downregulated/Wnt-associated cluster, vermillion for the upregulated/
    osteogenesis-associated cluster), scaled by the same global norm so
    darker = more significant within its own color family."""
    norm = (_sqrt[tf] - _lo) / (_hi - _lo)
    base = COL_WNT_HUB if tf in DOWN_GROUP else COL_OSTEO_HUB
    # blend from ~80%-white (least significant) to the full base color
    # (most significant) -- simple linear RGB interpolation
    r, g, b = int(base[1:3], 16), int(base[3:5], 16), int(base[5:7], 16)
    white_fraction = 0.15 + 0.65 * (1 - norm)  # 15% white (most significant) to 80% white (least significant)
    r2 = round(r + (255 - r) * white_fraction)
    g2 = round(g + (255 - g) * white_fraction)
    b2 = round(b + (255 - b) * white_fraction)
    return f"rgb({r2},{g2},{b2})"


# ---------------------------------------------------------------
# Node positions -- identical layout to tf_wnt_osteogenesis_network.py
# ---------------------------------------------------------------
nodes = {
    "Canonical<br>Wnt signaling": (0.0, 0.0, "hub_wnt"),
    "Osteogenesis":             (10.0, 0.0, "hub_osteo"),
    "LEF1":    (-2.3, 0.0, "tf"),
    "SP7":     (12.2, 2.4, "tf"),
    "RING1B":  (12.2, 0.8, "tf"),
    "NUCKS1":  (12.2, -0.8, "tf"),
    "IRF8":    (12.2, -2.4, "tf"),
    "SUZ12":   (5.0, 4.2, "tf"),
    "EZH2":    (5.0, 3.0, "tf"),
    "JARID2":  (5.0, 1.8, "tf"),
    "MTF2":    (5.0, 0.6, "tf"),
    "EP300":   (5.0, -0.6, "tf"),
    "BMI1":    (5.0, -1.8, "tf"),
    "SOX2":    (5.0, -3.0, "tf"),
    "PPARG":   (5.0, -4.2, "tf"),
}

# ---------------------------------------------------------------
# Edges -- identical literature basis to tf_wnt_osteogenesis_network.py
# ---------------------------------------------------------------
edges = [
    ("LEF1", "Canonical<br>Wnt signaling", "strong"),
    ("SUZ12", "Canonical<br>Wnt signaling", "indirect"),
    ("SUZ12", "Osteogenesis", "indirect"),
    ("EZH2", "Canonical<br>Wnt signaling", "indirect"),
    ("EZH2", "Osteogenesis", "indirect"),
    ("JARID2", "Canonical<br>Wnt signaling", "indirect"),
    ("JARID2", "Osteogenesis", "indirect"),
    ("MTF2", "Canonical<br>Wnt signaling", "indirect"),
    ("MTF2", "Osteogenesis", "indirect"),
    ("BMI1", "Canonical<br>Wnt signaling", "indirect"),
    ("BMI1", "Osteogenesis", "indirect"),
    ("RING1B", "Osteogenesis", "indirect"),
    ("SOX2", "Canonical<br>Wnt signaling", "indirect"),
    ("SOX2", "Osteogenesis", "antagonist"),
    ("SP7", "Osteogenesis", "strong"),
    ("EP300", "Canonical<br>Wnt signaling", "strong"),
    ("EP300", "Osteogenesis", "strong"),
    ("IRF8", "Osteogenesis", "weak"),
    ("PPARG", "Canonical<br>Wnt signaling", "antagonist"),
    ("PPARG", "Osteogenesis", "antagonist"),
    ("NUCKS1", "Osteogenesis", "weak"),
]

STYLE_MAP = {
    "strong":    dict(color=COL_EDGE_STRONG, dash="solid", width=2.0),
    "indirect":  dict(color=COL_EDGE_INDIRECT, dash="dash", width=1.4),
    "weak":      dict(color=COL_EDGE_WEAK, dash="dot", width=1.2),
    "antagonist": dict(color=COL_EDGE_ANTAGONIST, dash="solid", width=2.0),
}

fig = go.Figure()

for style_key in ["indirect", "weak", "strong", "antagonist"]:
    style = STYLE_MAP[style_key]
    xs, ys = [], []
    for src, dst, s in edges:
        if s != style_key:
            continue
        x0, y0, _ = nodes[src]
        x1, y1, _ = nodes[dst]
        xs += [x0, x1, None]
        ys += [y0, y1, None]
    label = {"strong": "Direct / well-established", "indirect": "Indirect (complex membership / context-dependent)",
              "weak": "Weak / associative (tentative)", "antagonist": "Antagonistic relationship"}[style_key]
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="lines",
        line=dict(color=style["color"], dash=style["dash"], width=style["width"]),
        name=label, hoverinfo="skip", showlegend=True,
    ))

for name, (x, y, kind) in nodes.items():
    if kind not in ("hub_wnt", "hub_osteo"):
        continue
    color = COL_WNT_HUB if kind == "hub_wnt" else COL_OSTEO_HUB
    fig.add_trace(go.Scatter(
        x=[x], y=[y], mode="markers+text", text=[name], textposition="middle center",
        textfont=dict(size=13, color="white", family="Liberation Sans"),
        marker=dict(size=95, color=color, line=dict(width=0)),
        hoverinfo="skip", showlegend=False,
    ))

tf_names = [n for n, (x, y, kind) in nodes.items() if kind == "tf"]
fig.add_trace(go.Scatter(
    x=[nodes[n][0] for n in tf_names],
    y=[nodes[n][1] for n in tf_names],
    mode="markers+text",
    text=tf_names,
    textposition="top center",
    textfont=dict(size=11, color="#1a1a1a", family="Liberation Sans"),
    marker=dict(
        size=[node_size(n) for n in tf_names],
        color=[node_color(n) for n in tf_names],
        line=dict(width=1.5, color=COL_NODE_BORDER),
    ),
    customdata=[[n, FDR[n]] for n in tf_names],
    hovertemplate="<b>%{customdata[0]}</b><br>FDR=%{customdata[1]:.2e}<extra></extra>",
    showlegend=False,
))

fig.update_layout(
    title=dict(text="Literature-supported TF\u2013pathway relationships overlaid with target-set enrichment evidence",
               font=dict(size=18, family="Liberation Sans", color="#1a1a1a")),
    plot_bgcolor="white", paper_bgcolor="white",
    xaxis=dict(visible=False, range=[-4, 14]),
    yaxis=dict(visible=False, range=[-5.5, 5.5], scaleanchor="x", scaleratio=1),
    legend=dict(orientation="h", yanchor="bottom", y=-0.07, xanchor="center", x=0.5, font=dict(size=10)),
    margin=dict(l=20, r=20, t=60, b=220),
    width=1050, height=920,
)

fig.add_annotation(
    text="Nodes = TF target sets from exploratory enrichment analysis. Node size/color = -log\u2081\u2080(BH-adjusted<br>"
         "enrichment p value), sqrt-scaled for visualization. Edges = literature-supported relationships with<br>"
         "canonical Wnt signaling or osteogenesis, independent of the enrichment statistics; edge style denotes<br>"
         "the curated evidence category (see tf_edge_evidence_table.csv). These associations do not establish<br>"
         "TF activity, direct binding, regulatory direction, or causality. No TF target set reached significance<br>"
         "in both queried databases.",
    xref="paper", yref="paper", x=0.5, y=-0.18, showarrow=False,
    align="center",
    font=dict(size=10, color="#666666", family="Liberation Sans"),
)

fig.write_html("tf_hybrid_network_fdr.html")
fig.write_image("tf_hybrid_network_fdr.png", scale=3)
print("Saved tf_hybrid_network_fdr.html and .png")
for n in tf_names:
    print(f"  {n:<8} FDR={FDR[n]:.2e}  size={node_size(n):.1f}  color={node_color(n)}")
