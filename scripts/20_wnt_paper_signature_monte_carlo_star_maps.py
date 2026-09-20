#!/usr/bin/env python3
"""Monte Carlo-supported Wnt/beta-catenin paper-signature star maps.

This script uses the available DESeq2 *summary* table and the four gene
signatures from Doumpas et al. (EMBO J 2019, doi:10.15252/embj.201798873).
It does NOT claim sample-level bootstrap or phenotype-label permutation,
because the summary table has no per-sample integer counts.

For each paper signature and current-study direction, it:
  1. selects the current top-2,000 genes by nominal p-value;
  2. calculates the observed signature overlap;
  3. draws baseMean-decile-matched random gene pools (Monte Carlo null);
  4. calculates an empirical p-value and BH FDR across all tests;
  5. scores each overlap gene across four sensitivity pools
     (p<0.05, p<0.10, top-2,000, top-5,000); and
  6. creates radial consensus star maps whose spoke width/opacity encodes
     threshold stability.

The human HEK293T paper signatures are contextual reference sets; overlap in
this mouse experiment is not independent validation of the paper mechanism.

Dependencies:
  python -m pip install pandas numpy matplotlib openpyxl

Example (repository root):
  python scripts/20_wnt_paper_signature_monte_carlo_star_maps.py \
    --deseq2 data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv \
    --signatures data/reference/Wnt_beta_catenin_LEFTCF_signatures.xlsx
"""

from __future__ import annotations

import argparse
import json
import math
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


PAPER_REFERENCE = "Doumpas et al., EMBO J 2019, doi:10.15252/embj.201798873"
N_MONTE_CARLO_DEFAULT = 10_000
RANDOM_SEED_DEFAULT = 20260917
TOP_N = 2_000
MAX_GENES_PER_MAP = 24

SIGNATURE_SHEETS = {
    "beta-catenin independent": "bCat indep (65 out of 231)",
    "beta-catenin dependent": "bCat dep (166 out of 231)",
    "TCF/LEF independent": "TCF indep (27 out of 166)",
    "TCF/LEF dependent": "TCF dep (139 out of 166)",
}

COLORS = {
    "downregulated": {"hub": "#0072B2", "gene": "#B9DDF0", "edge": "#4E8BA8"},
    "upregulated": {"hub": "#D55E00", "gene": "#F4C9A8", "edge": "#B87849"},
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--deseq2", type=Path, required=True, help="DESeq2 summary CSV/TXT")
    p.add_argument("--signatures", type=Path, required=True, help="Doumpas Table EV1 workbook")
    p.add_argument(
        "--outdir", type=Path,
        default=Path("results/wnt_paper_signature_monte_carlo"),
    )
    p.add_argument("--n-sim", type=int, default=N_MONTE_CARLO_DEFAULT)
    p.add_argument("--seed", type=int, default=RANDOM_SEED_DEFAULT)
    return p.parse_args()


def sniff_sep(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        line = handle.readline()
    return "\t" if line.count("\t") > line.count(",") else ","


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    x = p[ok]
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0, 1)
    out[np.where(ok)[0]] = adjusted
    return out


def load_deseq2(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep=sniff_sep(path))
    required = {"symbol", "baseMean", "log2FoldChange", "pvalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DESeq2 table missing columns: {sorted(missing)}")
    for column in ["baseMean", "log2FoldChange", "pvalue"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["symbol"] = df["symbol"].astype("string").str.strip().str.upper()
    df = df.dropna(subset=["symbol", "baseMean", "log2FoldChange", "pvalue"])
    df = (
        df.sort_values(["pvalue", "symbol"], kind="mergesort")
        .drop_duplicates("symbol", keep="first")
        .reset_index(drop=True)
    )
    if len(df) < TOP_N:
        raise ValueError(f"Only {len(df)} valid unique genes; need at least {TOP_N}.")
    df["padj_recomputed"] = bh_adjust(df["pvalue"])
    # Stable quantile bins, with duplicates handled explicitly.
    ranks = df["baseMean"].rank(method="first")
    df["expression_bin"] = pd.qcut(ranks, q=10, labels=False).astype(int)
    return df


def load_signatures(path: Path) -> dict[str, set[str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    available = pd.ExcelFile(path).sheet_names
    missing = set(SIGNATURE_SHEETS.values()) - set(available)
    if missing:
        raise ValueError(f"Signature workbook missing sheets: {sorted(missing)}")
    signatures: dict[str, set[str]] = {}
    for label, sheet in SIGNATURE_SHEETS.items():
        table = pd.read_excel(path, sheet_name=sheet)
        if table.shape[1] < 2:
            raise ValueError(f"Sheet {sheet!r} has no symbol column at position 2.")
        genes = set(
            table.iloc[:, 1].dropna().astype(str).str.strip().str.upper()
        )
        genes.discard("")
        signatures[label] = genes
    return signatures


def build_pools(df: pd.DataFrame) -> dict[str, dict[str, set[str]]]:
    scenarios = {
        "p005": df[df["pvalue"] < 0.05],
        "p010": df[df["pvalue"] < 0.10],
        "top2000": df.nsmallest(2_000, "pvalue"),
        "top5000": df.nsmallest(5_000, "pvalue"),
    }
    pools: dict[str, dict[str, set[str]]] = {}
    for scenario, sub in scenarios.items():
        pools[scenario] = {
            "downregulated": set(sub.loc[sub["log2FoldChange"] < 0, "symbol"]),
            "upregulated": set(sub.loc[sub["log2FoldChange"] > 0, "symbol"]),
        }
    observed = {
        d: len(pools["top2000"][d]) for d in ["downregulated", "upregulated"]
    }
    if observed != {"downregulated": 1428, "upregulated": 572}:
        raise RuntimeError(
            "Top-2,000 direction counts do not match the validated analysis: "
            f"{observed}"
        )
    return pools


def monte_carlo_test(
    universe: pd.DataFrame,
    observed_pool: set[str],
    signature: set[str],
    n_sim: int,
    rng: np.random.Generator,
) -> tuple[int, float, float, float, float]:
    observed = len(observed_pool & signature)
    # A random draw without replacement within each expression decile has a
    # hypergeometric overlap distribution. Drawing it directly is identical to
    # explicitly sampling thousands of gene-name sets, but is much faster and
    # avoids repeatedly filtering a large data frame.
    observed_rows = universe[universe["symbol"].isin(observed_pool)]
    null = np.zeros(n_sim, dtype=int)
    for expression_bin, draw_size in observed_rows["expression_bin"].value_counts().items():
        bin_rows = universe[universe["expression_bin"] == expression_bin]
        n_total = len(bin_rows)
        n_signature = int(bin_rows["symbol"].isin(signature).sum())
        null += rng.hypergeometric(
            ngood=n_signature,
            nbad=n_total - n_signature,
            nsample=int(draw_size),
            size=n_sim,
        )
    empirical_p = (1 + int(np.sum(null >= observed))) / (n_sim + 1)
    return (
        observed,
        empirical_p,
        float(null.mean()),
        float(np.quantile(null, 0.025)),
        float(np.quantile(null, 0.975)),
    )


def threshold_support(gene: str, direction: str, pools: dict) -> int:
    return sum(gene in pools[s][direction] for s in ["p005", "p010", "top2000", "top5000"])


def safe_name(value: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in value).strip("_")


def draw_star(
    row: pd.Series,
    genes: list[str],
    supports: dict[str, int],
    out_png: Path,
) -> plt.Figure:
    direction = row["direction"]
    color = COLORS[direction]
    shown = genes[:MAX_GENES_PER_MAP]
    n = len(shown)
    angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n, endpoint=False)
    radius = 1.0
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    for gene, angle in zip(shown, angles):
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        support = supports[gene]
        lw = 0.7 + 0.65 * support
        alpha = 0.25 + 0.18 * support
        ax.plot([0, x], [0, y], color=color["edge"], lw=lw, alpha=alpha, zorder=1)
        ax.text(
            x, y, gene, ha="center", va="center", fontsize=11.5,
            bbox=dict(
                boxstyle="round,pad=0.38,rounding_size=0.5",
                facecolor=color["gene"], edgecolor="#30343A", linewidth=1.3,
            ), zorder=3,
        )
    hub_edge = "#1B7837" if row["empirical_fdr_bh"] < 0.05 else "#30343A"
    ax.text(
        0, 0, row["signature"], ha="center", va="center", color="white",
        fontsize=15, fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.65,rounding_size=0.8",
            facecolor=color["hub"], edgecolor=hub_edge, linewidth=3.0,
        ), zorder=4,
    )
    title = f"{row['signature']} - {direction} consensus star map"
    subtitle = (
        f"Top-2,000 pool | observed overlap = {row['observed_overlap']} | "
        f"matched Monte Carlo p = {row['empirical_p']:.3g} | "
        f"BH FDR = {row['empirical_fdr_bh']:.3g}"
    )
    fig.suptitle(title, fontsize=22, fontweight="bold", color="#142B52", y=0.965)
    ax.set_title(subtitle, fontsize=13.5, color="#3D4A5C", pad=14)
    note = (
        "Spoke width and opacity encode inclusion across four gene-selection "
        "scenarios (p<0.05, p<0.10, top-2,000 and top-5,000). A green hub "
        "border indicates Monte Carlo BH FDR <0.05. Human HEK293T paper "
        "signatures are contextual and do not independently validate this "
        "mouse experiment."
    )
    fig.text(0.5, 0.035, textwrap.fill(note, 145), ha="center", va="bottom",
             fontsize=10.5, color="#555555")
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.23, 1.23)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0.06, right=0.94, top=0.84, bottom=0.13)
    fig.savefig(out_png, dpi=300, facecolor="white")
    return fig


def main() -> None:
    args = parse_args()
    if args.n_sim < 1000:
        raise ValueError("Use at least 1,000 Monte Carlo simulations.")
    args.outdir.mkdir(parents=True, exist_ok=True)
    map_dir = args.outdir / "individual_maps"
    map_dir.mkdir(exist_ok=True)

    df = load_deseq2(args.deseq2)
    signatures = load_signatures(args.signatures)
    pools = build_pools(df)
    rng = np.random.default_rng(args.seed)

    records = []
    for signature_name, signature_genes in signatures.items():
        detected_signature = signature_genes & set(df["symbol"])
        for direction in ["downregulated", "upregulated"]:
            pool = pools["top2000"][direction]
            observed, p_emp, null_mean, null_low, null_high = monte_carlo_test(
                df, pool, detected_signature, args.n_sim, rng
            )
            overlap = sorted(pool & detected_signature)
            records.append({
                "signature": signature_name,
                "direction": direction,
                "paper_signature_size": len(signature_genes),
                "detected_signature_size": len(detected_signature),
                "input_pool_size": len(pool),
                "observed_overlap": observed,
                "overlap_genes": ";".join(overlap),
                "empirical_p": p_emp,
                "null_mean": null_mean,
                "null_2.5pct": null_low,
                "null_97.5pct": null_high,
            })

    results = pd.DataFrame(records)
    results["empirical_fdr_bh"] = bh_adjust(results["empirical_p"])
    results = results.sort_values(["empirical_fdr_bh", "empirical_p", "signature", "direction"])
    results.to_csv(args.outdir / "wnt_signature_monte_carlo_results.csv", index=False)

    stability_rows = []
    figures = []
    for _, row in results.iterrows():
        genes = [g for g in str(row["overlap_genes"]).split(";") if g]
        if not genes:
            continue
        supports = {g: threshold_support(g, row["direction"], pools) for g in genes}
        genes = sorted(genes, key=lambda g: (-supports[g], g))
        for gene in genes:
            source = df.loc[df["symbol"].eq(gene)].iloc[0]
            stability_rows.append({
                "signature": row["signature"], "direction": row["direction"],
                "gene": gene, "threshold_support_count": supports[gene],
                "threshold_support_fraction": supports[gene] / 4,
                "log2FoldChange": source["log2FoldChange"],
                "pvalue": source["pvalue"],
                "padj_recomputed": source["padj_recomputed"],
            })
        stem = safe_name(f"{row['direction']}_{row['signature']}_monte_carlo_star")
        png = map_dir / f"{stem}.png"
        fig = draw_star(row, genes, supports, png)
        figures.append(fig)

    pd.DataFrame(stability_rows).to_csv(
        args.outdir / "wnt_signature_gene_stability.csv", index=False
    )
    combined_pdf = args.outdir / "wnt_paper_signature_monte_carlo_star_maps.pdf"
    with PdfPages(combined_pdf) as pdf:
        for fig in figures:
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "deseq2": str(args.deseq2), "signatures": str(args.signatures),
        "paper_reference": PAPER_REFERENCE, "n_sim": args.n_sim,
        "seed": args.seed, "matching": "baseMean deciles",
        "limitations": [
            "No sample-level counts were available.",
            "No sample bootstrap, leave-one-out DESeq2, or phenotype permutation was performed.",
            "Paper signatures are human HEK293T reference sets applied contextually to mouse data.",
        ],
    }
    (args.outdir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(f"Validated unique genes: {len(df)}")
    print("Top-2,000: 1428 downregulated, 572 upregulated")
    print(f"Monte Carlo tests: {len(results)} x {args.n_sim:,} simulations")
    print(f"Generated star maps: {len(figures)}")
    print(f"Combined PDF: {combined_pdf}")


if __name__ == "__main__":
    main()
