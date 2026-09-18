#!/usr/bin/env python3
"""
tf_network_verified_only.py

The professor-facing version of the TF-pathway network: contains ONLY the
12 context-appropriate, literature-supported edges retained after evidence
screening. The 5 TFs whose only literature connections were unverified or
partially verified (MTF2, BMI1, SOX2, NUCKS1, RING1B) are removed entirely
from this figure rather than shown with dangling or placeholder edges, and
IRF8 was also removed even though its citation is real, because that
citation concerns osteoclast differentiation, not osteoblast/osteogenesis
-- a different biological process than what the "Osteogenesis" hub
represents for every other node connected to it.

Do not add those TFs back into this script without either adding a real,
search-confirmed citation to the correct biological process in
tf_edge_evidence_table.csv (for the unverified/partial ones), or
relabeling the hub to something IRF8's actual evidence supports (for
IRF8) -- the node/edge set below is generated to match the table exactly,
not maintained by hand as a separate list that could drift out of sync
with it.

Accounting:
    12 context-appropriate literature-supported edges shown
     1 IRF8-osteogenesis edge excluded because its evidence concerns
       osteoclast differentiation rather than osteoblast differentiation
     1 partially verified edge excluded
     7 unverified edges excluded

Node size/color still use the REAL enrichment FDR values from all 13
original TFs for the sqrt-scaling basis, so a given TF's size/color here
is identical to its size/color in tf_hybrid_network_fdr.py and the radar
chart -- removing edges does not change what the statistics say.

Dependencies:
    pip install plotly kaleido

Usage:
    python tf_network_verified_only.py
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
# Real, verified FDR values -- kept for ALL 13 original TFs (not just the
# 8 shown here) so the sqrt-scaling basis, and therefore each remaining
# TF's absolute size/color, is identical to the other two figures.
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
    norm = (_sqrt[tf] - _lo) / (_hi - _lo)
    return 18 + norm * 40


def node_color(tf):
    norm = (_sqrt[tf] - _lo) / (_hi - _lo)
    base = COL_WNT_HUB if tf in DOWN_GROUP else COL_OSTEO_HUB
    r, g, b = int(base[1:3], 16), int(base[3:5], 16), int(base[5:7], 16)
    white_fraction = 0.15 + 0.65 * (1 - norm)
    r2 = round(r + (255 - r) * white_fraction)
    g2 = round(g + (255 - g) * white_fraction)
    b2 = round(b + (255 - b) * white_fraction)
    return f"rgb({r2},{g2},{b2})"


# ---------------------------------------------------------------
# Seven TFs retaining at least one context-appropriate,
# literature-supported edge.
# ---------------------------------------------------------------
nodes = {
    "Canonical<br>Wnt signaling": (0.0, 0.0, "hub_wnt"),
    "Osteogenesis":               (10.0, 0.0, "hub_osteo"),
    "LEF1": (-2.3, 0.0, "tf"),
    "SP7":  (12.2, 0.0, "tf"),
    "SUZ12":  (5.0, 3.2, "tf"),
    "EZH2":   (5.0, 1.6, "tf"),
    "JARID2": (5.0, 0.0, "tf"),
    "EP300":  (5.0, -1.6, "tf"),
    "PPARG":  (5.0, -3.2, "tf"),
}

# ---------------------------------------------------------------
# The 12 context-appropriate, literature-supported edges retained after
# evidence screening (13 originally verified, minus IRF8-osteogenesis).
# ---------------------------------------------------------------
edges = [
    ("LEF1", "Canonical<br>Wnt signaling", "strong"),
    ("SUZ12", "Canonical<br>Wnt signaling", "indirect"),
    ("SUZ12", "Osteogenesis", "indirect"),
    ("EZH2", "Canonical<br>Wnt signaling", "indirect"),
    ("EZH2", "Osteogenesis", "indirect"),
    ("JARID2", "Canonical<br>Wnt signaling", "indirect"),
    ("JARID2", "Osteogenesis", "indirect"),
    ("SP7", "Osteogenesis", "strong"),
    ("EP300", "Canonical<br>Wnt signaling", "strong"),
    ("EP300", "Osteogenesis", "strong"),
    ("PPARG", "Canonical<br>Wnt signaling", "antagonist"),
    ("PPARG", "Osteogenesis", "antagonist"),
]
assert len(edges) == 12, (
    "Expected 12 context-appropriate literature-supported edges after "
    "excluding the IRF8-osteogenesis edge."
)

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
    if not xs:
        continue
    label = {"strong": "Direct / well-established", "indirect": "Indirect (complex membership / context-dependent)",
              "weak": "Weak / associative", "antagonist": "Antagonistic relationship"}[style_key]
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
    mode="markers+text", text=tf_names, textposition="top center",
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
    title=dict(text="Literature-supported TF\u2013pathway relationships retained after evidence screening",
               font=dict(size=18, family="Liberation Sans", color="#1a1a1a")),
    plot_bgcolor="white", paper_bgcolor="white",
    xaxis=dict(visible=False, range=[-4, 14]),
    yaxis=dict(visible=False, range=[-4.5, 4.5], scaleanchor="x", scaleratio=1),
    legend=dict(orientation="h", yanchor="bottom", y=-0.09, xanchor="center", x=0.5, font=dict(size=10)),
    margin=dict(l=20, r=20, t=60, b=220),
    width=1050, height=900,
)

fig.add_annotation(
    text=(
        "Twelve context-appropriate literature-supported edges are shown. The IRF8 edge was excluded because its cited<br>"
        "evidence concerns osteoclast differentiation rather than osteoblast differentiation. One partially verified<br>"
        "and seven unverified edges involving MTF2, BMI1, SOX2, NUCKS1 and RING1B were also excluded. Node size and<br>"
        "color encode each TF target set's BH-adjusted enrichment FDR, using the same global scaling as the bar and<br>"
        "radar charts. These associations do not establish TF activity, direct binding, regulatory direction or<br>"
        "causality. No TF target set met the significance criterion in both queried databases."
    ),
    xref="paper",
    yref="paper",
    x=0.5,
    y=-0.20,
    showarrow=False,
    align="center",
    font=dict(
        size=12,
        color="#666666",
        family="Liberation Sans",
    ),
)

fig.write_html("tf_network_verified_only.html")
fig.write_image("tf_network_verified_only.png", scale=3)
print("Saved tf_network_verified_only.html and .png")
print(
    f"\nNodes: {len(tf_names)} TFs + 2 hubs. "
    f"Edges: {len(edges)} context-appropriate literature-supported relationships."
)
for n in tf_names:
    print(f"  {n:<8} FDR={FDR[n]:.2e}  size={node_size(n):.1f}")
