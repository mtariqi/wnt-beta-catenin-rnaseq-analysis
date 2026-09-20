#!/usr/bin/env python3
"""Create an orthology-verified, publication-style Doumpas signature figure.

The published human signatures are converted to mouse symbols only through the
consensus one-to-one orthology table.  A mouse ortholog is "represented" when
its symbol occurs in the valid DESeq2 differential-expression analysis
universe.  Representation is a coverage measure, not evidence of enrichment.

Expected final flows (published -> mapped -> represented):
  TCF/LEF dependent        139 -> 101 -> 97
  beta-catenin dependent   166 -> 120 -> 116
  beta-catenin independent  65 ->  43 -> 42
  TCF/LEF independent       27 ->  19 -> 19
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import openpyxl
import pandas as pd

PUBLICATION_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}
plt.rcParams.update(PUBLICATION_STYLE)


COLORS = {
    "NS": "#D3D3D3",
    "0.05 ≤ p < 0.10": "#D55E00",
    "p < 0.05": "#0072B2",
}
SIG_ORDER = ["NS", "0.05 ≤ p < 0.10", "p < 0.05"]
DIRECTION_MARKERS = {"Up": "o", "Down": "v"}
GUIDE = "#9E9E9E"
WNT_REFERENCE_GENES = {"AXIN2", "LEF1", "SP5", "NOTUM", "NKD1"}
WNT_REFERENCE_COLOR = "#E69F00"
N_EXPLORATORY_LABELS = 5

# Deterministic, panel-specific label placement in offset points.  These
# offsets are deliberately version-controlled instead of being optimized at
# render time, so SVG/PDF/PNG outputs reproduce identically across machines.
LABEL_OFFSETS = {
    "tcf_dep": {
        "MYH11": (7, 8),
        "ANKRD24": (-40, 6),
        "BCL11B": (10, -12),
        "HRH1": (-27, -8),
    },
    "bcat_dep": {
        "MYH11": (7, 9),
        "CHAC1": (8, 12),
    },
    "bcat_indep": {
        "COL6A5": (-38, -3),
        "PDE3B": (10, 9),
        "PRPH": (10, -8),
    },
}


PANELS = [
    {
        "key": "tcf_dep",
        "sheet": "TCF dep (139 out of 166)",
        "label": "TCF/LEF dependent",
        "published": 139,
        "mapped": 101,
        "represented": 97,
    },
    {
        "key": "bcat_dep",
        "sheet": "bCat dep (166 out of 231)",
        "label": "β-catenin dependent",
        "published": 166,
        "mapped": 120,
        "represented": 116,
    },
    {
        "key": "bcat_indep",
        "sheet": "bCat indep (65 out of 231)",
        "label": "β-catenin independent",
        "published": 65,
        "mapped": 43,
        "represented": 42,
    },
    {
        "key": "tcf_indep",
        "sheet": "TCF indep (27 out of 166)",
        "label": "TCF/LEF independent",
        "published": 27,
        "mapped": 19,
        "represented": 19,
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--deseq2",
        default="data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv",
    )
    parser.add_argument(
        "--signatures",
        default="data/reference/Wnt_beta_catenin_LEFTCF_signatures.xlsx",
    )
    parser.add_argument(
        "--ortholog-map",
        default="data/reference/human_mouse_consensus_one2one.csv",
    )
    parser.add_argument(
        "--outdir",
        default="results/wnt_signature_ranked_validation_verified_orthologs/nature_style",
    )
    parser.add_argument(
        "--basename",
        default="doumpas_signatures_orthology_verified_2x2",
    )
    parser.add_argument(
        "--allow-count-mismatch",
        action="store_true",
        help="Plot despite disagreement with the finalized audit counts.",
    )
    return parser.parse_args()


def require_file(path: str | Path) -> Path:
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"Required input not found: {resolved}")
    return resolved


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def classify_sig(pvalue: float) -> str:
    if pvalue < 0.05:
        return "p < 0.05"
    if pvalue < 0.10:
        return "0.05 ≤ p < 0.10"
    return "NS"


def load_deseq(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"symbol", "log2FoldChange", "pvalue"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"DESeq2 file is missing columns: {sorted(missing)}")

    df = df.copy()
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")

    # Inferential step: calculate fallback BH values across every row with a
    # valid gene-level p-value before removing records that cannot be displayed
    # by gene symbol. If the original DESeq2 padj column is available, preserve
    # it and use the fallback only where the source value is missing.
    valid_test = df["pvalue"].notna() & df["pvalue"].between(0, 1)
    df["padj_plot"] = np.nan
    if "padj" in df.columns:
        source_padj = pd.to_numeric(df["padj"], errors="coerce")
    else:
        source_padj = pd.Series(np.nan, index=df.index, dtype=float)
    recalculated = pd.Series(np.nan, index=df.index, dtype=float)
    recalculated.loc[valid_test] = bh_adjust(df.loc[valid_test, "pvalue"].to_numpy())
    df["padj_plot"] = source_padj.where(source_padj.notna(), recalculated)

    # Display step: require a symbol, fold change, and p-value, then select one
    # record per displayed symbol. This does not alter the BH testing universe.
    df = df.dropna(subset=["symbol", "log2FoldChange", "pvalue"])
    df = df.loc[(df["pvalue"] >= 0) & (df["pvalue"] <= 1)].copy()
    df["symbol_upper"] = df["symbol"].astype(str).str.strip().str.upper()
    df = df.loc[df["symbol_upper"].ne("")].copy()
    df = df.sort_values("pvalue").drop_duplicates("symbol_upper", keep="first")
    return df


def load_ortholog_map(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path)
    required = {"human_symbol", "mouse_symbol"}
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(f"Orthology map is missing columns: {sorted(missing)}")

    mapping = mapping.dropna(subset=["human_symbol", "mouse_symbol"]).copy()
    mapping["human_upper"] = mapping["human_symbol"].astype(str).str.strip().str.upper()
    mapping["mouse_upper"] = mapping["mouse_symbol"].astype(str).str.strip().str.upper()
    mapping = mapping.loc[
        mapping["human_upper"].ne("") & mapping["mouse_upper"].ne("")
    ].drop_duplicates(["human_upper", "mouse_upper"])

    duplicated_human = mapping["human_upper"].duplicated(keep=False)
    duplicated_mouse = mapping["mouse_upper"].duplicated(keep=False)
    if duplicated_human.any() or duplicated_mouse.any():
        raise ValueError(
            "The supplied orthology map is not one-to-one after normalization."
        )
    return mapping


def load_signature_genes(path: Path, sheet_name: str) -> list[str]:
    # Some source workbooks contain unsupported visual extensions. They do not
    # affect cell values used here.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Unknown extension is not supported and will be removed"
        )
        workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
        if sheet_name not in workbook.sheetnames:
            workbook.close()
            raise ValueError(f"Worksheet not found: {sheet_name}")
        worksheet = workbook[sheet_name]
        genes = []
        for row in worksheet.iter_rows(min_row=2, values_only=True):
            symbol = row[1] if len(row) > 1 else None
            if isinstance(symbol, str) and symbol.strip():
                genes.append(symbol.strip().upper())
        workbook.close()
    return list(dict.fromkeys(genes))


def prepare_signature(
    deseq: pd.DataFrame,
    orthologs: pd.DataFrame,
    human_genes: list[str],
) -> tuple[pd.DataFrame, int]:
    human_set = set(human_genes)
    mapped = orthologs.loc[orthologs["human_upper"].isin(human_set)].copy()
    mapped_mouse = set(mapped["mouse_upper"])
    sub = deseq.loc[deseq["symbol_upper"].isin(mapped_mouse)].copy()

    # Zero p-values can arise through output rounding; cap only for plotting.
    smallest_positive = np.nextafter(0.0, 1.0)
    sub["neglog10p"] = -np.log10(sub["pvalue"].clip(lower=smallest_positive))
    sub["Significance"] = sub["pvalue"].map(classify_sig)
    sub["Direction"] = np.where(sub["log2FoldChange"] > 0, "Up", "Down")
    sub["label_score"] = (
        sub["neglog10p"].pow(0.7) * sub["log2FoldChange"].abs()
    )
    return sub, len(mapped_mouse)


def validate_counts(
    panel: dict,
    human_genes: list[str],
    mapped_count: int,
    represented_count: int,
    allow_mismatch: bool,
) -> None:
    observed = (len(human_genes), mapped_count, represented_count)
    expected = (panel["published"], panel["mapped"], panel["represented"])
    if observed != expected:
        message = (
            f"{panel['label']} count mismatch: observed published/mapped/represented "
            f"{observed}, expected {expected}."
        )
        if allow_mismatch:
            warnings.warn(message, stacklevel=2)
        else:
            raise ValueError(message + " Use --allow-count-mismatch only for diagnostics.")


def select_labels(sub: pd.DataFrame) -> pd.DataFrame:
    # Always retain prespecified reference genes, irrespective of nominal
    # significance, then add only a small number of informative nominal hits.
    candidates = sub.loc[sub["pvalue"] < 0.10].copy()
    reference = sub.loc[sub["symbol_upper"].isin(WNT_REFERENCE_GENES)]
    remaining = candidates.loc[
        ~candidates["symbol_upper"].isin(WNT_REFERENCE_GENES)
    ].sort_values("label_score", ascending=False)
    return pd.concat(
        [reference, remaining.head(N_EXPLORATORY_LABELS)],
        ignore_index=False,
    ).drop_duplicates("symbol_upper")


def draw_panel(
    ax: plt.Axes,
    sub: pd.DataFrame,
    panel: dict,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    for sig_class in SIG_ORDER:
        for direction, marker in DIRECTION_MARKERS.items():
            mask = (sub["Significance"] == sig_class) & (
                sub["Direction"] == direction
            )
            if not mask.any():
                continue
            ax.scatter(
                sub.loc[mask, "log2FoldChange"],
                sub.loc[mask, "neglog10p"],
                s=34,
                marker=marker,
                facecolor=COLORS[sig_class],
                edgecolor="white",
                linewidth=0.35,
                alpha=0.88,
                zorder=3,
            )

    reference = sub["symbol_upper"].isin(WNT_REFERENCE_GENES)
    if reference.any():
        ax.scatter(
            sub.loc[reference, "log2FoldChange"],
            sub.loc[reference, "neglog10p"],
            s=82,
            facecolors="none",
            edgecolors=WNT_REFERENCE_COLOR,
            linewidths=1.15,
            zorder=6,
        )

    fdr_sig = sub["padj_plot"] < 0.05
    if fdr_sig.any():
        ax.scatter(
            sub.loc[fdr_sig, "log2FoldChange"],
            sub.loc[fdr_sig, "neglog10p"],
            s=72,
            marker="s",
            facecolors="none",
            edgecolors="black",
            linewidths=1.0,
            zorder=7,
        )

    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color=GUIDE, zorder=1)
    ax.axhline(-np.log10(0.10), ls=":", lw=0.7, color=GUIDE, zorder=1)
    ax.axvline(0, ls="--", lw=0.7, color=GUIDE, zorder=1)

    for _, row in select_labels(sub).iterrows():
        symbol = str(row["symbol_upper"])
        dx, dy = LABEL_OFFSETS.get(panel["key"], {}).get(symbol, (5, 4))
        displaced = (dx, dy) != (5, 4)
        ax.annotate(
            str(row["symbol"]),
            xy=(row["log2FoldChange"], row["neglog10p"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.2,
            fontstyle="italic",
            ha="left" if dx >= 0 else "right",
            va="bottom" if dy >= 0 else "top",
            arrowprops=(
                {"arrowstyle": "-", "lw": 0.4, "color": "#777777"}
                if displaced
                else None
            ),
            annotation_clip=True,
            zorder=8,
        )

    n05 = int((sub["pvalue"] < 0.05).sum())
    n05_to_10 = int(((sub["pvalue"] >= 0.05) & (sub["pvalue"] < 0.10)).sum())
    ax.set_title(
        f"{panel['label']}\n"
        f"{panel['represented']} represented | {panel['mapped']} mapped | "
        f"{panel['published']} published\n"
        f"nominal p<0.05: {n05} | 0.05≤p<0.10: {n05_to_10}",
        fontsize=9,
        fontweight="bold",
        pad=6,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r"log$_2$ fold change (mutant vs control)")
    ax.set_ylabel(r"$-log_{10}$(nominal p-value)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)


def build_legend(fig: plt.Figure, has_fdr_hits: bool) -> None:
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=COLORS[label],
            markeredgecolor="white",
            markersize=6,
            label=label,
        )
        for label in SIG_ORDER
    ]
    handles.extend(
        Line2D(
            [0],
            [0],
            marker=marker,
            linestyle="",
            markerfacecolor="#555555",
            markeredgecolor="white",
            markersize=6,
            label=direction,
        )
        for direction, marker in DIRECTION_MARKERS.items()
    )
    handles.extend(
        [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor="none",
                markeredgecolor=WNT_REFERENCE_COLOR,
                markersize=7,
                label="Prespecified canonical Wnt reference gene",
            ),
        ]
    )
    if has_fdr_hits:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="",
                markerfacecolor="none",
                markeredgecolor="black",
                markersize=6,
                label="BH FDR < 0.05",
            )
        )
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.072),
        ncol=4,
        frameon=False,
        handletextpad=0.4,
        columnspacing=1.2,
    )


def main() -> None:
    args = parse_args()
    deseq_path = require_file(args.deseq2)
    signatures_path = require_file(args.signatures)
    ortholog_path = require_file(args.ortholog_map)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    deseq = load_deseq(deseq_path)
    orthologs = load_ortholog_map(ortholog_path)
    prepared: dict[str, pd.DataFrame] = {}

    for panel in PANELS:
        human_genes = load_signature_genes(signatures_path, panel["sheet"])
        sub, mapped_count = prepare_signature(deseq, orthologs, human_genes)
        validate_counts(
            panel,
            human_genes,
            mapped_count,
            len(sub),
            args.allow_count_mismatch,
        )
        prepared[panel["key"]] = sub
        print(
            f"{panel['label']}: {len(human_genes)} published -> "
            f"{mapped_count} mapped -> {len(sub)} represented"
        )

    max_fc = max(float(sub["log2FoldChange"].abs().max()) for sub in prepared.values())
    max_y = max(float(sub["neglog10p"].max()) for sub in prepared.values())
    # Compact symmetric limits without clipping an observed gene.  A fixed
    # +/-3.8 range would hide the strongest negative fold changes here.
    x_bound = max(3.8, float(np.ceil((max_fc + 0.20) * 10.0) / 10.0))
    y_bound = max(1.0, float(np.ceil(max_y * 1.08)))
    xlim = (-x_bound, x_bound)
    ylim = (0.0, y_bound)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.5), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.085,
        right=0.985,
        # Reserve a distinct title band above the three-line panel headings.
        top=0.82,
        bottom=0.17,
        wspace=0.16,
        hspace=0.40,
    )

    for letter, ax, panel in zip("ABCD", axes.flat, PANELS):
        draw_panel(ax, prepared[panel["key"]], panel, xlim, ylim)
        ax.text(
            -0.14,
            1.04,
            letter,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
        )

    fig.suptitle(
        "Orthology-verified Wnt/β-catenin signatures",
        fontsize=11,
        fontweight="bold",
        y=0.975,
    )
    fig.text(
        0.5,
        0.942,
        "Current mutant-versus-control RNA-seq differential-expression analysis "
        "(n=4 control, n=4 mutant)",
        ha="center",
        va="top",
        fontsize=8.5,
        color="#333333",
    )
    has_fdr_hits = any(
        bool((sub["padj_plot"] < 0.05).any()) for sub in prepared.values()
    )
    build_legend(fig, has_fdr_hits)
    fig.text(
        0.5,
        0.018,
        "Mapped mouse orthologs were considered represented when present in the valid DESeq2 analysis universe. "
        "Signature sets overlap and counts are not additive. High representation indicates coverage, not pathway "
        "enrichment; orthology-verified ranked tests were non-significant after BH correction (minimum FDR=0.491).",
        ha="center",
        va="bottom",
        fontsize=7,
        color="#444444",
        wrap=True,
    )

    base = outdir / args.basename
    fig.savefig(base.with_suffix(".svg"), facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.savefig(base.with_suffix(".png"), dpi=600, facecolor="white")
    plt.close(fig)
    print(f"Created {base}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
