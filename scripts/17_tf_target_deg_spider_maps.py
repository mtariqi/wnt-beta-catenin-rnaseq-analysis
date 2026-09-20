#!/usr/bin/env python3
"""Create auditable TF-centred target-gene spider maps for Top-2,000 inputs.

The centre node is an enriched TF target set. Peripheral nodes are the exact
overlap genes reported in the Enrichr ``Genes`` field for that record. These
edges represent target-set membership, not TF activity, direct binding in the
studied tissue, regulatory direction, or causality.

Inputs (run from the repository root):
  results/tf_sensitivity/enrichment_full_by_scenario.csv
  results/tf_sensitivity/gene_pools/top2000_downregulated.csv
  results/tf_sensitivity/gene_pools/top2000_upregulated.csv

Outputs:
  results/tf_target_spider_maps/top2000_tf_target_spider_maps.pdf
  results/tf_target_spider_maps/top2000_tf_target_evidence_table.csv
  results/tf_target_spider_maps/top2000_tf_spider_selection_manifest.csv
  results/tf_target_spider_maps/individual_maps/*.png and *.svg

Dependencies: pandas, numpy, matplotlib, networkx
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
import numpy as np
import pandas as pd


ENRICHMENT = Path("results/tf_sensitivity/enrichment_full_by_scenario.csv")
POOL_DIR = Path("results/tf_sensitivity/gene_pools")
OUTDIR = Path("results/tf_target_spider_maps")
MAPDIR = OUTDIR / "individual_maps"

SCENARIO = "top2000"
FDR_CUTOFF = 0.05
TOP_N = {"downregulated": 8, "upregulated": 5}
EXPECTED_POOL_COUNTS = {"downregulated": 1428, "upregulated": 572}
MAX_GENES_PER_MAP = 15

BLUE = "#0072B2"
ORANGE = "#D55E00"
LIGHT_BLUE = "#B9DDF0"
LIGHT_ORANGE = "#F4C9A8"
DARK = "#222222"
GREY = "#737373"
WEB_GREY = "#D5D9DE"

# Canonical display label for the RNF2/RING1B synonym pair.
REGULATOR_ALIASES = {"RNF2": "RING1B"}

REQUIRED_COLUMNS = {
    "scenario", "direction", "Gene_set", "Term", "Adjusted P-value", "Genes"
}


def regulator_from_term(term: str) -> str:
    """Extract the leading regulator label while preserving the original term."""
    match = re.match(r"^\s*([A-Za-z0-9.-]+)", str(term))
    if not match:
        raise ValueError(f"Cannot extract regulator from term: {term!r}")
    return match.group(1).upper()


def split_overlap_genes(value: object) -> list[str]:
    """Parse Enrichr overlap genes; accept semicolon or comma separators."""
    if pd.isna(value):
        return []
    genes = [g.strip().upper() for g in re.split(r"[;,]", str(value))]
    return list(dict.fromkeys(g for g in genes if g))


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def load_pool(direction: str) -> set[str]:
    path = POOL_DIR / f"{SCENARIO}_{direction}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing validated gene pool: {path}")
    table = pd.read_csv(path)
    if "symbol" not in table.columns:
        raise ValueError(f"{path} must contain a 'symbol' column")
    genes = set(table["symbol"].dropna().astype(str).str.strip().str.upper())
    expected = EXPECTED_POOL_COUNTS[direction]
    if len(genes) != expected:
        raise RuntimeError(
            f"{direction} pool contains {len(genes)} unique genes; expected {expected}."
        )
    return genes


def load_and_select() -> tuple[pd.DataFrame, dict[str, set[str]]]:
    if not ENRICHMENT.exists():
        raise FileNotFoundError(
            f"Missing {ENRICHMENT}. Run scripts/04_threshold_deg_tf_sensitivity.py first."
        )

    full = pd.read_csv(ENRICHMENT)
    missing = REQUIRED_COLUMNS.difference(full.columns)
    if missing:
        raise ValueError(f"Enrichment table is missing columns: {sorted(missing)}")

    full["Adjusted P-value"] = pd.to_numeric(full["Adjusted P-value"], errors="coerce")
    full = full.loc[
        full["scenario"].eq(SCENARIO)
        & full["direction"].isin(TOP_N)
        & full["Adjusted P-value"].lt(FDR_CUTOFF)
    ].copy()
    if full.empty:
        raise RuntimeError("No significant Top-2,000 enrichment records were found.")

    full["regulator"] = full.get(
        "regulator", full["Term"].map(regulator_from_term)
    ).astype(str).str.upper().replace(REGULATOR_ALIASES)
    full["overlap_genes"] = full["Genes"].map(split_overlap_genes)
    full["overlap_gene_count"] = full["overlap_genes"].map(len)

    pools = {direction: load_pool(direction) for direction in TOP_N}
    for idx, row in full.iterrows():
        outside = sorted(set(row["overlap_genes"]) - pools[row["direction"]])
        if outside:
            raise RuntimeError(
                f"Overlap mismatch for {row['regulator']} ({row['direction']}): "
                f"{outside[:10]} are absent from its validated Top-2,000 pool."
            )

    # One record per regulator/direction across both libraries: retain the
    # smallest enrichment FDR, with deterministic tie breaking.
    ranked = full.sort_values(
        ["direction", "regulator", "Adjusted P-value", "Gene_set", "Term"],
        kind="mergesort",
    )
    best = ranked.drop_duplicates(["direction", "regulator"], keep="first")

    selected_parts = []
    for direction, n_records in TOP_N.items():
        part = (
            best.loc[best["direction"].eq(direction)]
            .sort_values(["Adjusted P-value", "regulator"], kind="mergesort")
            .head(n_records)
            .copy()
        )
        if len(part) != n_records:
            raise RuntimeError(
                f"Expected {n_records} unique {direction} TFs, found {len(part)}."
            )
        part["selection_rank"] = np.arange(1, len(part) + 1)
        selected_parts.append(part)

    selected = pd.concat(selected_parts, ignore_index=True)
    if selected.duplicated(["direction", "regulator"]).any():
        raise RuntimeError("Duplicate regulator-direction selections detected.")
    if len(selected) != sum(TOP_N.values()):
        raise RuntimeError("Expected exactly 13 selected spider-map records.")
    return selected, pools


def radial_positions(tf: str, genes: list[str]) -> dict[str, tuple[float, float]]:
    positions = {tf: (0.0, 0.0)}
    angles = np.linspace(0, 2 * math.pi, len(genes), endpoint=False) + math.pi / 2
    for gene, angle in zip(genes, angles):
        positions[gene] = (math.cos(angle), math.sin(angle))
    return positions


def draw_spider_web(ax: plt.Axes, genes: list[str]) -> None:
    """Draw polygonal guide rings; these are not biological gene-gene edges."""
    angles = np.linspace(0, 2 * math.pi, len(genes), endpoint=False) + math.pi / 2
    for radius, width, alpha in [(0.34, 1.0, 0.60), (0.67, 1.2, 0.75), (1.0, 1.5, 0.90)]:
        xs = [radius * math.cos(angle) for angle in angles]
        ys = [radius * math.sin(angle) for angle in angles]
        ax.plot(xs + [xs[0]], ys + [ys[0]], color=WEB_GREY,
                linewidth=width, alpha=alpha, zorder=0)
    for angle in angles:
        ax.plot([0, math.cos(angle)], [0, math.sin(angle)],
                color=WEB_GREY, linewidth=0.9, alpha=0.60, zorder=0)


def draw_map(row: pd.Series) -> plt.Figure:
    tf = row["regulator"]
    direction = row["direction"]
    all_genes = row["overlap_genes"]
    # A TF can occur in its own enrichment overlap. Keep that evidence
    # in the table, but do not let it overwrite the central TF position.
    peripheral_genes = [gene for gene in all_genes if gene != tf]
    genes = peripheral_genes[:MAX_GENES_PER_MAP]
    if not genes:
        raise RuntimeError(f"{tf} has no parsed overlap genes.")

    hub_color = BLUE if direction == "downregulated" else ORANGE
    gene_color = LIGHT_BLUE if direction == "downregulated" else LIGHT_ORANGE
    graph = nx.Graph()
    graph.add_node(tf, node_type="TF")
    graph.add_edges_from((tf, gene) for gene in genes)
    pos = radial_positions(tf, genes)
    hub_size = 4200 + max(0, len(tf) - 4) * 1100

    fig, ax = plt.subplots(figsize=(13.33, 7.5), constrained_layout=False)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    draw_spider_web(ax, genes)
    nx.draw_networkx_edges(graph, pos, ax=ax, width=2.4, edge_color=GREY, alpha=0.95)
    nx.draw_networkx_nodes(
        graph, pos, nodelist=genes, ax=ax, node_size=1800,
        node_color=gene_color, edgecolors="#555555", linewidths=1.8,
    )
    nx.draw_networkx_nodes(
        graph, pos, nodelist=[tf], ax=ax, node_size=hub_size,
        node_color=hub_color, edgecolors="#333333", linewidths=2.2,
    )
    nx.draw_networkx_labels(
        graph, pos, labels={g: g for g in genes}, ax=ax,
        font_size=15, font_weight="bold", font_color=DARK,
    )
    nx.draw_networkx_labels(
        graph, pos, labels={tf: tf}, ax=ax,
        font_size=18, font_weight="bold", font_color="white",
    )

    direction_label = direction.replace("regulated", "regulated candidate genes")
    fig.suptitle(
        f"{tf}-centered target-gene spider map - Top-2,000 {direction_label}",
        fontsize=20, fontweight="bold", y=0.965,
    )
    ax.set_title(
        f"{row['Gene_set']}  |  enrichment FDR = {row['Adjusted P-value']:.2e}  |  "
        f"overlap genes = {len(all_genes)}",
        fontsize=16, color="#374151", pad=14,
    )
    if len(peripheral_genes) > len(genes):
        self_note = " (self-overlap retained in table)" if tf in all_genes else ""
        display_note = (f"Showing {len(genes)} of {len(peripheral_genes)} peripheral overlap genes"
                        f"{self_note}; complete list is in the evidence table.")
    else:
        display_note = f"Showing all {len(genes)} exact overlap genes."
    fig.text(0.5, 0.055, display_note, ha="center", fontsize=13, color="#4B5563")
    fig.text(
        0.5, 0.025,
        "Spokes indicate Enrichr overlap; pale web rings are layout guides only. "
        "No TF activity, direct binding, gene-to-gene link, or causality is established.",
        ha="center", fontsize=11.5, color="#4B5563",
    )
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-1.22, 1.22)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.84, bottom=0.11)
    return fig


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    MAPDIR.mkdir(parents=True, exist_ok=True)
    selected, _ = load_and_select()

    manifest_rows = []
    evidence_rows = []
    pdf_path = OUTDIR / "top2000_tf_target_spider_maps.pdf"

    with PdfPages(pdf_path) as pdf:
        for _, row in selected.sort_values(["direction", "selection_rank"]).iterrows():
            tf = row["regulator"]
            direction = row["direction"]
            library = row["Gene_set"]
            stem = safe_name(f"{direction}_{tf}_{library}")
            png_path = MAPDIR / f"{stem}.png"
            svg_path = MAPDIR / f"{stem}.svg"

            fig = draw_map(row)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
            fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
            plt.close(fig)

            shown = [gene for gene in row["overlap_genes"] if gene != tf][:MAX_GENES_PER_MAP]
            manifest_rows.append({
                "scenario": SCENARIO,
                "direction": direction,
                "selection_rank": int(row["selection_rank"]),
                "regulator": tf,
                "Gene_set": library,
                "Term": row["Term"],
                "Adjusted P-value": row["Adjusted P-value"],
                "overlap_gene_count": len(row["overlap_genes"]),
                "displayed_gene_count": len(shown),
                "displayed_genes": ";".join(shown),
                "map_generated": True,
                "png_file": str(png_path),
                "svg_file": str(svg_path),
            })
            for gene in row["overlap_genes"]:
                evidence_rows.append({
                    "scenario": SCENARIO,
                    "direction": direction,
                    "regulator": tf,
                    "Gene_set": library,
                    "Term": row["Term"],
                    "Adjusted P-value": row["Adjusted P-value"],
                    "overlap_gene": gene,
                    "selected_for_figure": gene in shown,
                })

    manifest = pd.DataFrame(manifest_rows)
    evidence = pd.DataFrame(evidence_rows)
    manifest_path = OUTDIR / "top2000_tf_spider_selection_manifest.csv"
    evidence_path = OUTDIR / "top2000_tf_target_evidence_table.csv"
    manifest.to_csv(manifest_path, index=False)
    evidence.to_csv(evidence_path, index=False)

    assert len(manifest) == 13
    assert manifest["map_generated"].all()
    assert manifest.groupby("direction").size().to_dict() == TOP_N

    print("Top-2,000 selection verified: 2000 total (1428 down, 572 up)")
    print(f"Generated {len(manifest)} TF-centred maps (8 downregulated-list, 5 upregulated-list).")
    print(f"Saved combined PDF: {pdf_path}")
    print(f"Saved selection manifest: {manifest_path}")
    print(f"Saved complete TF-overlap evidence table: {evidence_path}")


if __name__ == "__main__":
    main()
