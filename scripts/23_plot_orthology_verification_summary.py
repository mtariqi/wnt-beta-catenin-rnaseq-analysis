#!/usr/bin/env python3
"""Create a publication-ready summary of verified human-to-mouse orthology.

The script reads the finalized outputs from
``scripts/22_verify_human_mouse_orthology.py`` and writes PNG, PDF, and SVG
figures. Signature sets overlap, so per-signature counts are intentionally
presented separately from the overall unique-gene mapping outcome.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orthology-dir",
        type=Path,
        default=Path("results/wnt_orthology_verification"),
        help="Directory containing the finalized orthology CSV/JSON outputs.",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path(
            "results/wnt_orthology_verification/orthology_verification_summary"
        ),
        help="Output path without an extension.",
    )
    return parser.parse_args()


def load_verified_results(orthology_dir: Path) -> tuple[pd.DataFrame, dict]:
    retention_path = orthology_dir / "orthology_retention_by_signature.csv"
    provenance_path = orthology_dir / "orthology_provenance.json"

    if not retention_path.exists():
        raise FileNotFoundError(retention_path)
    if not provenance_path.exists():
        raise FileNotFoundError(provenance_path)

    retention = pd.read_csv(retention_path)
    required = {
        "signature",
        "paper_genes",
        "accepted_consensus_one2one",
        "retention_fraction",
    }
    missing = required.difference(retention.columns)
    if missing:
        raise ValueError(f"Missing retention columns: {sorted(missing)}")

    with provenance_path.open(encoding="utf-8") as handle:
        provenance = json.load(handle)

    return retention, provenance


def make_figure(retention: pd.DataFrame, provenance: dict, out_prefix: Path) -> None:
    evaluated = int(provenance.get("unique_human_signature_genes", 231))
    accepted = int(provenance["accepted_unique_pairs"])
    excluded = evaluated - accepted

    if evaluated <= 0 or accepted < 0 or accepted > evaluated:
        raise ValueError("Invalid overall orthology counts in provenance file")

    preferred_order = [
        "TCF/LEF dependent",
        "TCF/LEF independent",
        "beta-catenin dependent",
        "beta-catenin independent",
    ]
    order = {name: index for index, name in enumerate(preferred_order)}
    retention = retention.copy()
    retention["_order"] = retention["signature"].map(order).fillna(999)
    retention = retention.sort_values("_order").reset_index(drop=True)
    retention["retention_pct"] = retention["retention_fraction"] * 100

    display_names = (
        retention["signature"]
        .str.replace("beta-catenin", "β-catenin", regex=False)
        .tolist()
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
        }
    )
    fig = plt.figure(figsize=(14, 7.5), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.45], wspace=0.34)
    ax_overall = fig.add_subplot(grid[0, 0])
    ax_signature = fig.add_subplot(grid[0, 1])

    navy = "#18345F"
    accepted_color = "#2E86AB"
    excluded_color = "#D9E2EC"
    bar_color = "#4B8F8C"
    edge_color = "#263746"

    wedges, _ = ax_overall.pie(
        [accepted, excluded],
        startangle=90,
        counterclock=False,
        colors=[accepted_color, excluded_color],
        wedgeprops={"width": 0.34, "edgecolor": "white", "linewidth": 2},
    )
    ax_overall.text(
        0,
        0.07,
        f"{accepted}/{evaluated}",
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        color=navy,
    )
    ax_overall.text(
        0,
        -0.17,
        f"{accepted / evaluated:.1%} retained",
        ha="center",
        va="center",
        fontsize=12,
        color="#44546A",
    )
    ax_overall.set_title("Overall unique-gene mapping", color=navy, pad=16)
    ax_overall.legend(
        wedges,
        [f"Accepted consensus one-to-one: {accepted}", f"Excluded/unresolved: {excluded}"],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.19),
        frameon=False,
        fontsize=10,
    )

    y = range(len(retention))
    bars = ax_signature.barh(
        y,
        retention["retention_pct"],
        color=bar_color,
        edgecolor=edge_color,
        linewidth=0.8,
        height=0.62,
    )
    ax_signature.set_yticks(list(y), display_names)
    ax_signature.invert_yaxis()
    ax_signature.set_xlim(0, 100)
    ax_signature.set_xlabel("Consensus one-to-one ortholog retention (%)")
    ax_signature.set_title("Retention within published Wnt signatures", color=navy, pad=16)
    ax_signature.grid(axis="x", color="#D9E2EC", linewidth=0.8)
    ax_signature.set_axisbelow(True)
    ax_signature.spines[["top", "right", "left"]].set_visible(False)

    for bar, (_, row) in zip(bars, retention.iterrows()):
        accepted_n = int(row["accepted_consensus_one2one"])
        paper_n = int(row["paper_genes"])
        percentage = float(row["retention_pct"])
        ax_signature.text(
            min(percentage + 1.2, 93),
            bar.get_y() + bar.get_height() / 2,
            f"{accepted_n}/{paper_n} ({percentage:.1f}%)",
            va="center",
            ha="left",
            fontsize=10.5,
            fontweight="bold",
            color=edge_color,
        )

    fig.suptitle(
        "Human-to-mouse orthology verification of published Wnt gene signatures",
        fontsize=18,
        fontweight="bold",
        color=navy,
        y=0.98,
    )
    fig.subplots_adjust(left=0.07, right=0.97, top=0.86, bottom=0.18)
    fig.text(
        0.5,
        0.025,
        "Acceptance required a consensus one-to-one assignment under the Ensembl–MGI rule. "
        "Signature sets overlap; per-signature counts are not additive. PLAT and CALCOCO1 "
        "remained query-unresolved in the finalized run.",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#4A5568",
        wrap=True,
    )

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            out_prefix.with_suffix(f".{extension}"),
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    retention, provenance = load_verified_results(args.orthology_dir)
    make_figure(retention, provenance, args.out_prefix)
    print(f"Created {args.out_prefix}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
