#!/usr/bin/env python3
"""Plot canonical-axis and alternative-response Wnt signature trends.

This visualization deliberately avoids equating TCF/LEF-independent or
beta-catenin-independent transcription with a proven non-canonical Wnt
pathway. It reads the finalized orthology-verified Monte Carlo results and
labels the biology according to Doumpas et al. (2018).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(
            "results/wnt_signature_ranked_validation_verified_orthologs/"
            "ranked_signature_monte_carlo.csv"
        ),
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path(
            "results/wnt_signature_ranked_validation_verified_orthologs/"
            "canonical_axis_alternative_response_summary"
        ),
    )
    return parser.parse_args()


def load_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    data = pd.read_csv(path)
    required = {
        "signature",
        "detected_size",
        "nes",
        "empirical_p_two_sided",
        "fdr_bh_two_sided",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    data = data.copy()
    data["signature_key"] = (
        data["signature"].str.lower().str.replace("β", "beta", regex=False)
    )
    if data["signature_key"].duplicated().any():
        raise ValueError("Expected one row per signature")
    return data.set_index("signature_key")


def row_for(data: pd.DataFrame, signature: str) -> pd.Series:
    key = signature.lower().replace("β", "beta")
    if key not in data.index:
        raise ValueError(f"Signature not found: {signature}")
    return data.loc[key]


def plot_panel(ax, rows, title, subtitle, color):
    labels = [label for label, _ in rows]
    values = [float(row["nes"]) for _, row in rows]
    bars = ax.barh(
        range(len(rows)),
        values,
        color=color,
        edgecolor="#263746",
        linewidth=0.8,
        height=0.58,
    )
    ax.set_yticks(range(len(rows)), labels)
    ax.invert_yaxis()
    ax.axvline(0, color="#263746", linewidth=1.1)
    ax.set_xlim(-1.45, 0.15)
    ax.grid(axis="x", color="#D9E2EC", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title(title, loc="left", color="#18345F", fontsize=14, pad=27)
    ax.text(
        0,
        1.04,
        subtitle,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="#52616B",
    )
    for bar, (_, row) in zip(bars, rows):
        nes = float(row["nes"])
        p_value = float(row["empirical_p_two_sided"])
        fdr = float(row["fdr_bh_two_sided"])
        detected = int(row["detected_size"])
        ax.text(
            nes + 0.035,
            bar.get_y() + bar.get_height() / 2,
            f"NES={nes:.3f}\nP={p_value:.3f}; FDR={fdr:.3f}; n={detected}",
            ha="left",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color="white" if nes < -0.95 else "#263746",
        )


def make_figure(data: pd.DataFrame, out_prefix: Path) -> None:
    canonical_rows = [
        ("β-catenin dependent", row_for(data, "beta-catenin dependent")),
        ("TCF/LEF dependent", row_for(data, "TCF/LEF dependent")),
    ]
    alternative_rows = [
        ("TCF/LEF independent", row_for(data, "TCF/LEF independent")),
        ("β-catenin independent", row_for(data, "beta-catenin independent")),
    ]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.5), sharex=True)
    fig.patch.set_facecolor("white")

    plot_panel(
        axes[0],
        canonical_rows,
        "A. Canonical-axis dependency signatures",
        "Dependence on β-catenin or TCF/LEF supports the canonical transcriptional axis.",
        "#2E86AB",
    )
    plot_panel(
        axes[1],
        alternative_rows,
        "B. Alternative or independent response signatures",
        "Independence indicates an alternative response class, not proof of non-canonical Wnt.",
        "#D77A2B",
    )
    for ax in axes:
        ax.set_xlabel("Expression-matched normalized enrichment score (NES)")

    fig.suptitle(
        "Published Wnt signature trends in the current orthology-verified RNA-seq results",
        fontsize=18,
        fontweight="bold",
        color="#18345F",
        y=0.97,
    )
    fig.text(
        0.5,
        0.105,
        "All four signatures trend toward lower expression; none is statistically significant "
        "after BH correction (minimum FDR = 0.491).",
        ha="center",
        fontsize=11.5,
        fontweight="bold",
        color="#7A1F1F",
    )
    fig.text(
        0.5,
        0.045,
        "Interpretation follows Doumpas et al. (2018): TCF/LEF-independent β-catenin activity may "
        "represent an alternative branch within canonical Wnt/β-catenin signaling. These data do "
        "not establish activation or suppression of a ligand-defined non-canonical Wnt pathway.",
        ha="center",
        va="bottom",
        fontsize=9.4,
        color="#52616B",
        wrap=True,
    )
    fig.subplots_adjust(left=0.12, right=0.98, top=0.80, bottom=0.20, wspace=0.28)

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
    data = load_results(args.results)
    make_figure(data, args.out_prefix)
    print(f"Created {args.out_prefix}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
