#!/usr/bin/env python3
"""Create an audited TF threshold-stability heatmap and exact agreement tables."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SCENARIOS = ["p005", "p010", "top2000", "top5000"]
SCENARIO_LABELS = {
    "p005": "p < 0.05", "p010": "p < 0.10",
    "top2000": "Top 2,000", "top5000": "Top 5,000",
}
DATABASES = ["ChEA_2022", "TRRUST_Transcription_Factors_2019"]
DATABASE_LABELS = {
    "ChEA_2022": "ChEA", "TRRUST_Transcription_Factors_2019": "TRRUST",
}
PRIORITY = [
    "SUZ12", "MTF2", "EZH2", "SOX2", "JARID2", "RING1B", "LEF1", "BMI1",
    "SP7", "NUCKS1", "EP300", "IRF8", "PPARG", "RUNX2", "SMAD1",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=Path("results/tf_sensitivity/combined_all_enrichment_results.csv"))
    p.add_argument("--outdir", type=Path, default=Path("results/tf_sensitivity"))
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--top-extra", type=int, default=10)
    return p.parse_args()


def parse_regulator(term: str) -> str:
    match = re.match(r"^([^ _]+)", str(term))
    return match.group(1).upper() if match else str(term).upper()


def main() -> None:
    a = parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.input)
    required = {"scenario", "direction", "Gene_set", "Term", "Adjusted P-value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if "regulator" not in df:
        df["regulator"] = df["Term"].map(parse_regulator)
    df["regulator"] = df["regulator"].astype(str).str.upper()
    df = df[df["scenario"].isin(SCENARIOS) & df["Gene_set"].isin(DATABASES)].copy()

    # Minimum adjusted p across same-regulator terms, with the exact winning term retained.
    best = (
        df.sort_values("Adjusted P-value")
        .drop_duplicates(["scenario", "direction", "Gene_set", "regulator"])
        .copy()
    )
    best["FDR_significant"] = best["Adjusted P-value"] < a.fdr
    best["minus_log10_FDR"] = -np.log10(best["Adjusted P-value"].clip(lower=1e-300))

    # True dual-database agreement requires the same scenario, direction, and parsed regulator.
    agreement = (
        best[best["FDR_significant"]]
        .groupby(["scenario", "direction", "regulator"], as_index=False)
        .agg(
            significant_database_count=("Gene_set", "nunique"),
            significant_databases=("Gene_set", lambda x: ";".join(sorted(set(x)))),
        )
    )
    dual = agreement[agreement["significant_database_count"] == len(DATABASES)].copy()
    dual_details = best.merge(
        dual[["scenario", "direction", "regulator"]],
        on=["scenario", "direction", "regulator"], how="inner",
    )
    dual_details = dual_details[dual_details["FDR_significant"]].sort_values(
        ["scenario", "direction", "regulator", "Gene_set", "Adjusted P-value"]
    )
    dual.to_csv(a.outdir / "dual_database_regulators_by_scenario.csv", index=False)
    dual_details.to_csv(a.outdir / "dual_database_significant_details.csv", index=False)

    recurrence = (
        best.groupby(["direction", "Gene_set", "regulator"], as_index=False)
        .agg(
            scenarios_tested=("scenario", "nunique"),
            scenarios_FDR_significant=("FDR_significant", "sum"),
            minimum_adjusted_p=("Adjusted P-value", "min"),
        )
    )
    recurrence["stable_all_four"] = recurrence["scenarios_FDR_significant"].eq(4)
    recurrence.to_csv(a.outdir / "audited_regulator_threshold_stability.csv", index=False)

    extra = (
        recurrence.sort_values(["scenarios_FDR_significant", "minimum_adjusted_p"], ascending=[False, True])
        .loc[lambda x: ~x["regulator"].isin(PRIORITY), "regulator"]
        .drop_duplicates().head(a.top_extra).tolist()
    )
    regulators = [r for r in PRIORITY if r in set(best["regulator"])] + extra
    columns = [
        (s, d, db) for s in SCENARIOS
        for d in ["downregulated", "upregulated"] for db in DATABASES
    ]
    matrix = pd.DataFrame(index=regulators)
    for scenario, direction, database in columns:
        label = f"{SCENARIO_LABELS[scenario]}\n{direction[:4].title()} | {DATABASE_LABELS[database]}"
        series = best[
            best["scenario"].eq(scenario)
            & best["direction"].eq(direction)
            & best["Gene_set"].eq(database)
        ].set_index("regulator")["minus_log10_FDR"]
        matrix[label] = series.reindex(regulators)

    matrix.to_csv(a.outdir / "threshold_stability_heatmap_values.csv", index_label="regulator")
    annotations = matrix.map(
        lambda value: "●" if pd.notna(value) and value > -np.log10(a.fdr) else ""
    )
    sns.set_theme(style="white", font_scale=0.85)
    fig_width = max(15, 0.85 * matrix.shape[1])
    fig_height = max(8, 0.38 * matrix.shape[0] + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        matrix, cmap="Blues", vmin=0,
        mask=matrix.isna(), linewidths=0.35, linecolor="white",
        annot=annotations, fmt="", annot_kws={"fontsize": 8, "color": "black"},
        cbar_kws={"label": "−log₁₀(BH-adjusted enrichment p value)"}, ax=ax,
    )
    ax.set_title("TF target-set enrichment stability across input definitions", fontsize=17, weight="bold", pad=14)
    ax.set_xlabel("Input definition, gene-list direction and annotation database")
    ax.set_ylabel("Parsed target-set regulator")
    ax.tick_params(axis="x", labelrotation=55, labelsize=8)
    ax.tick_params(axis="y", labelrotation=0)
    fig.text(
        0.5, 0.01,
        "Filled circle: BH-adjusted enrichment p < 0.05. Gene-list thresholds use nominal DESeq2 p values; "
        "enrichment FDR is database-, direction- and scenario-specific.",
        ha="center", fontsize=9, color="#555555",
    )
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(a.outdir / "tf_threshold_stability_heatmap.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(a.outdir / "tf_threshold_stability_heatmap.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Dual-database scenario-direction-regulator hits: {len(dual)}")
    print(dual.to_string(index=False) if len(dual) else "None")
    print(f"\nSaved heatmap and audited tables to: {a.outdir}")


if __name__ == "__main__":
    main()
