#!/usr/bin/env python3
"""Professor-approved validation of published Wnt/beta-catenin signatures.

The analysis uses the four gene lists and published WT summary statistics in
Doumpas et al. (EMBO J, 2019; doi:10.15252/embj.201798873) as reference data.
It does not download or reanalyse the paper's raw sequencing reads/counts.

Outputs
-------
1. Symbol-compatibility audit between the human paper lists and mouse DESeq2.
2. Expression-decile-matched Monte Carlo enrichment of each signature in a
   signed, genome-wide DESeq2 ranking.
3. Threshold sensitivity across nominal-P and top-N gene pools, by direction.
4. Cross-study direction concordance using paper WT logFC versus current log2FC.
5. Publication-ready summary figures and a machine-readable provenance file.

Important interpretation
------------------------
Upper-case symbol matches are candidate human-mouse ortholog-compatible pairs,
not formal orthology proof. Supply --ortholog-map for a curated mapping table
with columns human_symbol and mouse_symbol when one becomes available.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, spearmanr


REFERENCE = "Doumpas et al. (2019), EMBO Journal 38:e98873"
DOI = "https://doi.org/10.15252/embj.201798873"
DEFAULT_SEED = 20260917
DEFAULT_N_SIM = 10_000

SHEETS = {
    "beta-catenin independent": "bCat indep (65 out of 231)",
    "beta-catenin dependent": "bCat dep (166 out of 231)",
    "TCF/LEF independent": "TCF indep (27 out of 166)",
    "TCF/LEF dependent": "TCF dep (139 out of 166)",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--deseq2", required=True, type=Path)
    p.add_argument("--signatures", required=True, type=Path)
    p.add_argument("--outdir", type=Path,
                   default=Path("results/wnt_signature_ranked_validation"))
    p.add_argument("--ortholog-map", type=Path, default=None,
                   help="Optional CSV/TSV with human_symbol and mouse_symbol")
    p.add_argument("--n-sim", type=int, default=DEFAULT_N_SIM)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return p.parse_args()


def sniff_sep(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        first = handle.readline()
    return "\t" if first.count("\t") > first.count(",") else ","


def bh_adjust(values: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(pd.to_numeric(values, errors="coerce"), dtype=float)
    out = np.full(p.size, np.nan)
    ok = np.isfinite(p)
    x = p[ok]
    if not x.size:
        return out
    order = np.argsort(x)
    adjusted_sorted = x[order] * x.size / np.arange(1, x.size + 1)
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = np.clip(adjusted_sorted, 0, 1)
    out[np.flatnonzero(ok)] = adjusted
    return out


def load_deseq2(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep=sniff_sep(path))
    required = {"symbol", "baseMean", "log2FoldChange", "pvalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DESeq2 file missing columns: {sorted(missing)}")
    for col in ["baseMean", "log2FoldChange", "pvalue"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["symbol"] = df["symbol"].astype("string").str.strip().str.upper()
    df = df.dropna(subset=["symbol", "baseMean", "log2FoldChange", "pvalue"])
    df = (df.sort_values(["pvalue", "symbol"], kind="mergesort")
            .drop_duplicates("symbol", keep="first").reset_index(drop=True))
    df["pvalue_clipped"] = df["pvalue"].clip(lower=np.nextafter(0, 1))
    # A signed significance ranking is transparent and reproducible when the
    # Wald statistic is unavailable from a summary export.
    df["rank_score"] = (np.sign(df["log2FoldChange"])
                        * -np.log10(df["pvalue_clipped"]).clip(upper=300))
    ranks = df["baseMean"].rank(method="first")
    df["expression_bin"] = pd.qcut(ranks, 10, labels=False).astype(int)
    return df


def load_ortholog_map(path: Path | None) -> dict[str, str] | None:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(path)
    tab = pd.read_csv(path, sep=sniff_sep(path))
    required = {"human_symbol", "mouse_symbol"}
    if not required.issubset(tab.columns):
        raise ValueError("Ortholog map needs human_symbol and mouse_symbol columns")
    tab = tab.dropna(subset=list(required)).copy()
    tab["human_symbol"] = tab["human_symbol"].astype(str).str.strip().str.upper()
    tab["mouse_symbol"] = tab["mouse_symbol"].astype(str).str.strip().str.upper()
    # One-to-many mappings require biological adjudication; keep only unique
    # human-to-mouse assignments to prevent silent duplication.
    counts = tab.groupby("human_symbol")["mouse_symbol"].nunique()
    unique_human = set(counts[counts.eq(1)].index)
    tab = tab[tab["human_symbol"].isin(unique_human)].drop_duplicates("human_symbol")
    return dict(zip(tab["human_symbol"], tab["mouse_symbol"]))


def load_signatures(path: Path, ortholog: dict[str, str] | None) -> dict[str, pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(path)
    available = set(pd.ExcelFile(path).sheet_names)
    missing = set(SHEETS.values()) - available
    if missing:
        raise ValueError(f"Workbook missing sheets: {sorted(missing)}")
    output = {}
    for label, sheet in SHEETS.items():
        tab = pd.read_excel(path, sheet_name=sheet)
        needed = {"symbol", "logFC WT", "pValue WT"}
        if not needed.issubset(tab.columns):
            raise ValueError(f"{sheet!r} missing columns: {sorted(needed-set(tab.columns))}")
        tab = tab[["symbol", "logFC WT", "pValue WT"]].copy()
        tab.columns = ["human_symbol", "paper_logfc_wt", "paper_pvalue_wt"]
        tab["human_symbol"] = tab["human_symbol"].astype("string").str.strip().str.upper()
        tab["paper_logfc_wt"] = pd.to_numeric(tab["paper_logfc_wt"], errors="coerce")
        tab["paper_pvalue_wt"] = pd.to_numeric(tab["paper_pvalue_wt"], errors="coerce")
        tab = tab.dropna(subset=["human_symbol"]).drop_duplicates("human_symbol")
        if ortholog is None:
            tab["mouse_symbol"] = tab["human_symbol"]
            tab["mapping_status"] = "symbol-compatible; orthology not formally verified"
        else:
            tab["mouse_symbol"] = tab["human_symbol"].map(ortholog)
            tab["mapping_status"] = np.where(tab["mouse_symbol"].notna(),
                                             "curated ortholog map", "unmapped")
        tab["signature"] = label
        output[label] = tab
    return output


def monte_carlo_rank_test(df: pd.DataFrame, signature: set[str], n_sim: int,
                          rng: np.random.Generator) -> dict[str, float]:
    members = df[df["symbol"].isin(signature)]
    if len(members) < 2:
        return dict(observed_mean=np.nan, null_mean=np.nan, null_sd=np.nan,
                    nes=np.nan, empirical_p_two_sided=np.nan,
                    empirical_p_up=np.nan, empirical_p_down=np.nan)
    observed = float(members["rank_score"].mean())
    null = np.zeros(n_sim, dtype=float)
    total_draw = 0
    for bin_id, draw_n in members["expression_bin"].value_counts().items():
        scores = df.loc[df["expression_bin"].eq(bin_id), "rank_score"].to_numpy()
        draw_n = int(draw_n)
        total_draw += draw_n
        # Independent samples within each simulation; no replacement within a
        # sampled set. This preserves the signature's expression-bin profile.
        for i in range(n_sim):
            null[i] += rng.choice(scores, size=draw_n, replace=False).sum()
    null /= total_draw
    center = float(null.mean())
    sd = float(null.std(ddof=1))
    nes = (observed - center) / sd if sd > 0 else np.nan
    two = (1 + int(np.sum(np.abs(null-center) >= abs(observed-center)))) / (n_sim+1)
    up = (1 + int(np.sum(null >= observed))) / (n_sim+1)
    down = (1 + int(np.sum(null <= observed))) / (n_sim+1)
    return dict(observed_mean=observed, null_mean=center, null_sd=sd, nes=nes,
                empirical_p_two_sided=two, empirical_p_up=up,
                empirical_p_down=down)


def ranked_analysis(df: pd.DataFrame, signatures: dict[str, pd.DataFrame],
                    n_sim: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    universe = set(df["symbol"])
    for label, tab in signatures.items():
        mapped = set(tab["mouse_symbol"].dropna())
        detected = mapped & universe
        result = monte_carlo_rank_test(df, detected, n_sim, rng)
        rows.append({"signature": label, "paper_size": len(tab),
                     "mapped_size": len(mapped), "detected_size": len(detected),
                     **result})
    out = pd.DataFrame(rows)
    out["fdr_bh_two_sided"] = bh_adjust(out["empirical_p_two_sided"])
    return out.sort_values(["fdr_bh_two_sided", "empirical_p_two_sided"])


def threshold_analysis(df: pd.DataFrame, signatures: dict[str, pd.DataFrame]) -> pd.DataFrame:
    scenarios = {
        "p<0.01": df[df["pvalue"] < 0.01],
        "p<0.05": df[df["pvalue"] < 0.05],
        "p<0.10": df[df["pvalue"] < 0.10],
        "top-1000": df.nsmallest(1000, "pvalue"),
        "top-2000": df.nsmallest(2000, "pvalue"),
        "top-5000": df.nsmallest(5000, "pvalue"),
    }
    universe = set(df["symbol"])
    rows = []
    for scenario, selected in scenarios.items():
        for direction, sign in [("downregulated", -1), ("upregulated", 1)]:
            pool = set(selected.loc[np.sign(selected["log2FoldChange"]).eq(sign), "symbol"])
            for label, tab in signatures.items():
                sig = set(tab["mouse_symbol"].dropna()) & universe
                a = len(pool & sig)
                b = len(pool - sig)
                c = len(sig - pool)
                d = len(universe - pool - sig)
                odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                rows.append({"scenario": scenario, "direction": direction,
                             "signature": label, "pool_size": len(pool),
                             "detected_signature_size": len(sig), "overlap": a,
                             "odds_ratio": odds, "fisher_p": p,
                             "overlap_genes": ";".join(sorted(pool & sig))})
    out = pd.DataFrame(rows)
    out["fdr_bh"] = bh_adjust(out["fisher_p"])
    return out.sort_values(["fdr_bh", "fisher_p"])


def mapping_and_concordance(df: pd.DataFrame,
                            signatures: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    current = df[["symbol", "log2FoldChange", "pvalue", "baseMean"]].rename(
        columns={"symbol": "mouse_symbol", "log2FoldChange": "current_log2fc",
                 "pvalue": "current_pvalue"})
    audits, summaries = [], []
    for label, tab in signatures.items():
        merged = tab.merge(current, on="mouse_symbol", how="left")
        merged["detected_current"] = merged["current_log2fc"].notna()
        merged["direction_concordant"] = np.where(
            merged["detected_current"] & merged["paper_logfc_wt"].notna(),
            np.sign(merged["paper_logfc_wt"]) == np.sign(merged["current_log2fc"]),
            np.nan)
        audits.append(merged)
        paired = merged.dropna(subset=["paper_logfc_wt", "current_log2fc"])
        n = len(paired)
        concordant = int(paired["direction_concordant"].astype(bool).sum()) if n else 0
        rho, rho_p = spearmanr(paired["paper_logfc_wt"], paired["current_log2fc"]) if n >= 3 else (np.nan, np.nan)
        binom_p = binomtest(concordant, n, 0.5, alternative="two-sided").pvalue if n else np.nan
        summaries.append({"signature": label, "paper_size": len(tab),
                          "detected_pairs": n, "direction_concordant": concordant,
                          "direction_concordance_fraction": concordant/n if n else np.nan,
                          "binomial_p_two_sided": binom_p,
                          "spearman_rho": rho, "spearman_p": rho_p})
    summary = pd.DataFrame(summaries)
    summary["binomial_fdr_bh"] = bh_adjust(summary["binomial_p_two_sided"])
    summary["spearman_fdr_bh"] = bh_adjust(summary["spearman_p"])
    return pd.concat(audits, ignore_index=True), summary


def make_figures(ranked: pd.DataFrame, thresholds: pd.DataFrame,
                 concordance: pd.DataFrame, outdir: Path) -> None:
    colors = ["#3B82A0" if x < 0 else "#D97732" for x in ranked["nes"]]
    fig, ax = plt.subplots(figsize=(10, 5.8))
    y = np.arange(len(ranked))
    ax.barh(y, ranked["nes"], color=colors, edgecolor="#30343A")
    ax.set_yticks(y, ranked["signature"])
    ax.axvline(0, color="#30343A", lw=1)
    ax.set_xlabel("Expression-matched normalized enrichment score")
    ax.set_title("Published Wnt signatures in the current ranked RNA-seq results",
                 weight="bold", color="#142B52")
    for i, (_, r) in enumerate(ranked.iterrows()):
        x = r["nes"]
        # Place the statistic inside the bar so long signature labels and the
        # annotation never collide in slide-sized exports.
        ax.text(x / 2, i, f"FDR={r['fdr_bh_two_sided']:.3f}",
                ha="center", va="center", fontsize=9.5, color="white",
                fontweight="bold")
    ax.set_xlim(min(-1.4, float(ranked["nes"].min()) - 0.12), 0.12)
    fig.subplots_adjust(left=0.34, right=0.97, top=0.88, bottom=0.14)
    fig.savefig(outdir / "ranked_signature_enrichment.png", dpi=300, facecolor="white")
    plt.close(fig)

    top = thresholds[thresholds["scenario"].eq("top-2000")].copy()
    top["label"] = top["signature"] + "\n" + top["direction"]
    fig, ax = plt.subplots(figsize=(12, 6.2))
    vals = np.log2(top["odds_ratio"].replace([np.inf, 0], np.nan))
    ax.bar(np.arange(len(top)), vals, color=["#3B82A0" if d == "downregulated" else "#D97732" for d in top["direction"]])
    ax.set_xticks(np.arange(len(top)), top["label"], rotation=35, ha="right")
    ax.axhline(0, color="#30343A", lw=1)
    ax.set_ylabel("log2 odds ratio")
    ax.set_title("Top-2,000 directional overlap effect sizes", weight="bold", color="#142B52")
    fig.tight_layout()
    fig.savefig(outdir / "top2000_overlap_effect_sizes.png", dpi=300, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.8))
    x = np.arange(len(concordance))
    ax.bar(x, concordance["direction_concordance_fraction"], color="#6A8E3A", edgecolor="#30343A")
    ax.axhline(0.5, color="#777777", ls="--", label="Chance expectation")
    ax.set_ylim(0, 1)
    ax.set_xticks(x, concordance["signature"], rotation=25, ha="right")
    ax.set_ylabel("Direction-concordant fraction")
    ax.set_title("Paper WT response versus current-study fold-change direction",
                 weight="bold", color="#142B52")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "cross_study_direction_concordance.png", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.n_sim < 1000:
        raise ValueError("Use at least 1,000 simulations")
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = load_deseq2(args.deseq2)
    ortholog = load_ortholog_map(args.ortholog_map)
    signatures = load_signatures(args.signatures, ortholog)

    audit, concordance = mapping_and_concordance(df, signatures)
    ranked = ranked_analysis(df, signatures, args.n_sim, args.seed)
    thresholds = threshold_analysis(df, signatures)

    audit.to_csv(args.outdir / "paper_signature_symbol_mapping_audit.csv", index=False)
    concordance.to_csv(args.outdir / "cross_study_direction_concordance.csv", index=False)
    ranked.to_csv(args.outdir / "ranked_signature_monte_carlo.csv", index=False)
    thresholds.to_csv(args.outdir / "threshold_sensitivity.csv", index=False)
    make_figures(ranked, thresholds, concordance, args.outdir)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference": REFERENCE, "doi": DOI,
        "paper_data_used": "four gene lists plus published WT logFC and pValue columns",
        "paper_raw_data_reanalysed": False,
        "deseq2": str(args.deseq2), "signatures": str(args.signatures),
        "ortholog_map": str(args.ortholog_map) if args.ortholog_map else None,
        "mapping_interpretation": ("curated map supplied" if ortholog else
            "uppercase symbol compatibility only; formal orthology not claimed"),
        "n_sim": args.n_sim, "seed": args.seed,
        "rank_metric": "sign(log2FoldChange) * -log10(nominal pvalue), capped at 300",
        "monte_carlo_matching": "baseMean deciles",
        "multiplicity": "Benjamini-Hochberg within each output family",
    }
    (args.outdir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8")

    print(f"Valid unique current-study genes: {len(df):,}")
    print(f"Ranked tests: {len(ranked)} signatures x {args.n_sim:,} simulations")
    print(f"Threshold tests: {len(thresholds)}")
    print(f"Output directory: {args.outdir}")
    print("Raw paper RNA-seq data reanalysed: no")


if __name__ == "__main__":
    main()
