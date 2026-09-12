#!/usr/bin/env python3
"""Sensitivity analysis for upstream TF target-set enrichment.

Compares four pre-specified gene-list definitions:
  * nominal p < 0.05
  * nominal p < 0.10
  * 2,000 smallest nominal p-values
  * 5,000 smallest nominal p-values

Gene selection uses nominal DESeq2 p-values without a fold-change magnitude
cutoff, as requested. Up- and downregulated genes are analyzed separately.
BH-adjusted p-values reported here are Enrichr term-level values within each
library/direction/threshold analysis; they are not DESeq2-adjusted gene p-values.

The Enrichr call requires internet access. Use --audit-only to generate and
validate all input lists without submitting them.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LIBRARIES = ["TRRUST_Transcription_Factors_2019", "ChEA_2022"]
COLORS = {"downregulated": "#0072B2", "upregulated": "#D55E00"}
REQUIRED = {"symbol", "pvalue", "log2FoldChange"}
SCENARIOS = [
    ("p005", "Nominal p < 0.05", "p", 0.05),
    ("p010", "Nominal p < 0.10", "p", 0.10),
    ("top2000", "Top 2,000 by nominal p", "top", 2000),
    ("top5000", "Top 5,000 by nominal p", "top", 5000),
]


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg adjustment with monotonicity and original order."""
    values = np.asarray(pvalues, dtype=float)
    if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must be finite and between 0 and 1.")
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    return restored


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--outdir", type=Path, default=Path("results/tf_sensitivity"))
    p.add_argument("--organism", default="mouse")
    p.add_argument("--enrichment-fdr", type=float, default=0.05)
    p.add_argument("--audit-only", action="store_true")
    p.add_argument(
        "--resume", action="store_true",
        help="Reuse a scenario's existing all-results CSV files when present.",
    )
    return p.parse_args()


def read_input(path: Path) -> tuple[pd.DataFrame, dict]:
    sep = "\t" if path.suffix.lower() in {".txt", ".tsv"} else ","
    raw = pd.read_csv(path, sep=sep)
    missing = REQUIRED - set(raw.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    valid = raw.dropna(subset=list(REQUIRED)).copy()
    valid["symbol"] = valid["symbol"].astype(str).str.strip().str.upper()
    valid = valid[valid["symbol"].ne("") & valid["log2FoldChange"].ne(0)].copy()
    # Deterministic symbol-level representative.
    valid = (
        valid.sort_values(["pvalue", "symbol"], kind="stable")
        .drop_duplicates("symbol", keep="first")
        .reset_index(drop=True)
    )
    valid["BH_FDR_AUDIT"] = bh_adjust(valid["pvalue"])
    audit = {
        "source_rows": int(len(raw)),
        "valid_unique_nonzero_direction_symbols": int(len(valid)),
        "DESeq2_BH_FDR_lt_0.05_symbols": int((valid["BH_FDR_AUDIT"] < 0.05).sum()),
        "important_note": (
            "BH_FDR_AUDIT is recomputed across available symbol-level nominal p-values "
            "and is not substituted for a DESeq2-native padj column."
        ),
    }
    return valid, audit


def select(table: pd.DataFrame, mode: str, value: float) -> pd.DataFrame:
    if mode == "p":
        chosen = table[table["pvalue"] < value].copy()
    elif mode == "top":
        n = int(value)
        if n > len(table):
            raise ValueError(f"Requested top {n:,}, but only {len(table):,} symbols are available.")
        chosen = table.nsmallest(n, "pvalue", keep="first").copy()
    else:
        raise ValueError(mode)
    chosen["direction"] = np.where(
        chosen["log2FoldChange"] > 0, "upregulated", "downregulated"
    )
    if set(chosen["direction"]) != set(COLORS):
        raise ValueError(f"{mode}={value} produced an empty direction.")
    return chosen


def regulator(term: str) -> str:
    match = re.match(r"^([^ _]+)", str(term))
    return match.group(1).upper() if match else str(term).upper()


def run_one(
    genes: list[str], scenario: str, direction: str, organism: str,
    outdir: Path, resume: bool,
) -> pd.DataFrame:
    import gseapy as gp

    scenario_dir = outdir / scenario
    scenario_dir.mkdir(parents=True, exist_ok=True)
    result_file = scenario_dir / f"{direction}_all_results.csv"
    if resume and result_file.exists():
        result = pd.read_csv(result_file)
    else:
        enr = gp.enrichr(
            gene_list=genes, gene_sets=LIBRARIES, organism=organism,
            outdir=str(scenario_dir / f"enrichr_{direction}"),
            cutoff=1.0, no_plot=True,
        )
        result = enr.results.copy()
        if result.empty:
            raise RuntimeError(f"No Enrichr results for {scenario}/{direction}")
        result.to_csv(result_file, index=False)
    needed = {"Gene_set", "Term", "P-value", "Adjusted P-value", "Odds Ratio"}
    missing = needed - set(result.columns)
    if missing:
        raise RuntimeError(f"Enrichr output missing columns: {sorted(missing)}")
    result["scenario"] = scenario
    result["direction"] = direction
    result["input_gene_count"] = len(genes)
    result["regulator"] = result["Term"].map(regulator)
    result["mouse_labelled"] = result["Term"].astype(str).str.contains(
        r"\bmouse\b", case=False, regex=True
    )
    return result


def plot_input_counts(counts: pd.DataFrame, output: Path) -> None:
    pivot = counts.pivot(index="scenario_label", columns="direction", values="gene_count")
    pivot = pivot.loc[[x[1] for x in SCENARIOS]]
    ax = pivot[["downregulated", "upregulated"]].plot.bar(
        figsize=(10, 6), color=[COLORS["downregulated"], COLORS["upregulated"]],
        edgecolor="black", linewidth=0.5,
    )
    ax.set_title("Input-gene counts across sensitivity scenarios", fontsize=17, weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Unique gene symbols")
    ax.legend(["Downregulated", "Upregulated"], frameon=False)
    ax.tick_params(axis="x", rotation=0)
    ax.spines[["top", "right"]].set_visible(False)
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=9, fmt="%d")
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def plot_significant_counts(summary: pd.DataFrame, output: Path) -> None:
    labels = {x[0]: x[1] for x in SCENARIOS}
    data = summary.copy()
    data["scenario_label"] = data["scenario"].map(labels)
    data["series"] = data["direction"].str.title() + " | " + data["Gene_set"].str.replace("_", " ")
    pivot = data.pivot(index="scenario_label", columns="series", values="significant_records").fillna(0)
    pivot = pivot.loc[[x[1] for x in SCENARIOS]]
    ax = pivot.plot.bar(figsize=(12, 6), width=0.82)
    ax.set_title("FDR-significant TF target-set records across input definitions", fontsize=16, weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Enrichment records with BH-adjusted p < 0.05")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def main() -> None:
    a = args()
    if not 0 < a.enrichment_fdr < 1:
        raise ValueError("--enrichment-fdr must be between 0 and 1.")
    a.outdir.mkdir(parents=True, exist_ok=True)
    table, audit = read_input(a.input)

    registry = {}
    gseapy_version = None
    if not a.audit_only:
        try:
            import gseapy as gp
        except ImportError as exc:
            raise RuntimeError(
                "GSEApy is required for enrichment. Install dependencies with "
                "'python -m pip install gseapy pandas numpy matplotlib'."
            ) from exc
        gseapy_version = gp.__version__
        available = set(gp.get_library_name(organism=a.organism))
        registry = {lib: lib in available for lib in LIBRARIES}
        if not all(registry.values()):
            raise RuntimeError(f"Required Enrichr library unavailable: {registry}")

    count_rows, combined = [], []
    for key, label, mode, value in SCENARIOS:
        chosen = select(table, mode, value)
        scenario_dir = a.outdir / key
        scenario_dir.mkdir(parents=True, exist_ok=True)
        chosen[["symbol", "log2FoldChange", "pvalue", "BH_FDR_AUDIT", "direction"]].to_csv(
            scenario_dir / "selected_gene_audit.csv", index=False
        )
        for direction in ["downregulated", "upregulated"]:
            genes = sorted(chosen.loc[chosen["direction"] == direction, "symbol"].unique())
            pd.Series(genes, name="symbol").to_csv(
                scenario_dir / f"{direction}_input_genes.csv", index=False
            )
            count_rows.append({
                "scenario": key, "scenario_label": label, "selection_mode": mode,
                "selection_value": value, "direction": direction,
                "gene_count": len(genes),
            })
            if not a.audit_only:
                combined.append(run_one(
                    genes, key, direction, a.organism, a.outdir, a.resume
                ))

    counts = pd.DataFrame(count_rows)
    counts.to_csv(a.outdir / "threshold_input_counts.csv", index=False)
    plot_input_counts(counts, a.outdir / "threshold_input_counts.png")

    metadata = {
        "python_version": sys.version.split()[0],
        "gseapy_version": gseapy_version,
        "organism": a.organism,
        "libraries": LIBRARIES,
        "library_registry_check": registry,
        "gene_selection_uses_nominal_p_without_abs_log2fc_cutoff": True,
        "enrichment_fdr_cutoff": a.enrichment_fdr,
        "audit_only": a.audit_only,
        **audit,
    }

    if not a.audit_only:
        all_results = pd.concat(combined, ignore_index=True)
        all_results.to_csv(a.outdir / "combined_all_enrichment_results.csv", index=False)
        sig = all_results[all_results["Adjusted P-value"] < a.enrichment_fdr].copy()
        sig.to_csv(a.outdir / "combined_FDR05_enrichment_results.csv", index=False)

        summary = (
            all_results.assign(significant=all_results["Adjusted P-value"] < a.enrichment_fdr)
            .groupby(["scenario", "direction", "Gene_set"], as_index=False)
            .agg(
                tested_records=("Term", "size"),
                significant_records=("significant", "sum"),
                minimum_adjusted_p=("Adjusted P-value", "min"),
            )
        )
        summary.to_csv(a.outdir / "scenario_enrichment_summary.csv", index=False)
        plot_significant_counts(summary, a.outdir / "scenario_significant_record_counts.png")

        # One row per regulator, direction, and library. Stability is defined as
        # the number of pre-specified scenarios in which any matching term is FDR-significant.
        stability = (
            all_results.assign(significant=all_results["Adjusted P-value"] < a.enrichment_fdr)
            .groupby(["direction", "Gene_set", "regulator"], as_index=False)
            .agg(
                scenarios_tested=("scenario", "nunique"),
                scenarios_FDR_significant=("significant", "sum"),
                minimum_adjusted_p=("Adjusted P-value", "min"),
                significant_scenarios=(
                    "scenario",
                    lambda x: ";".join(sorted(set(x[all_results.loc[x.index, "Adjusted P-value"] < a.enrichment_fdr]))),
                ),
            )
        )
        # Sum above may count multiple terms within a scenario; compute exact scenario recurrence.
        exact = (
            sig.groupby(["direction", "Gene_set", "regulator"])["scenario"]
            .nunique().rename("scenarios_FDR_significant_exact")
        )
        stability = stability.merge(exact, on=["direction", "Gene_set", "regulator"], how="left")
        stability["scenarios_FDR_significant"] = stability.pop(
            "scenarios_FDR_significant_exact"
        ).fillna(0).astype(int)
        stability["stable_all_four"] = stability["scenarios_FDR_significant"].eq(len(SCENARIOS))
        stability = stability.sort_values(
            ["stable_all_four", "scenarios_FDR_significant", "minimum_adjusted_p"],
            ascending=[False, False, True],
        )
        stability.to_csv(a.outdir / "regulator_threshold_stability.csv", index=False)

        # Cross-database replication at each scenario and direction.
        per_db = (
            sig.groupby(["scenario", "direction", "regulator"])["Gene_set"]
            .agg(lambda x: ";".join(sorted(set(x))))
            .reset_index(name="significant_libraries")
        )
        per_db["significant_in_both_databases"] = per_db["significant_libraries"].map(
            lambda x: set(x.split(";")).issuperset(LIBRARIES)
        )
        per_db.to_csv(a.outdir / "cross_database_agreement_by_scenario.csv", index=False)
        metadata["FDR_significant_records_total_across_scenarios"] = int(len(sig))
        metadata["regulator_library_direction_combinations_stable_all_four"] = int(
            stability["stable_all_four"].sum()
        )
        metadata["scenario_direction_regulators_significant_in_both_databases"] = int(
            per_db["significant_in_both_databases"].sum()
        )

    (a.outdir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print("\nINPUT COUNTS\n")
    print(counts.to_string(index=False))
    print(f"\nSaved sensitivity analysis to: {a.outdir}")
    if a.audit_only:
        print("Audit-only mode: no Enrichr requests were submitted.")


if __name__ == "__main__":
    main()
