#!/usr/bin/env python3
"""
04_threshold_deg_tf_sensitivity.py  (v2 -- corrected per review)

Consolidates Professor Hsu's first request: generate nominal-p-selected
candidate-gene pools at four
selection thresholds, split each by direction, run upstream TF-target
enrichment (ChEA_2022 + TRRUST_2019) for every pool/direction combination,
and compare TF stability across thresholds and databases.

This script consolidates the previously developed threshold-sensitivity
workflow into one reproducible implementation. (Not described as
"independently verified" beyond what is actually true: the candidate-gene-pool
generation logic was tested against the real DESeq2 table and matches
previously established counts exactly; the live Enrichr calls could not be
executed in the environment this script was written in and have not been
run end-to-end -- verify the enrichment half's real output yourself on
first run.)

PIPELINE:
    1. Load the real DESeq2 table (tab-delimited .txt by default -- see
       CANDIDATE_DESEQ2_PATHS; delimiter is auto-detected, not assumed).
    2. Deterministically deduplicate by symbol, keeping the row with the
       strongest nominal evidence (lowest p-value, then largest |log2FC|,
       then largest baseMean) rather than whichever row happens to appear
       first in the file.
    3. Generate four gene pools: nominal p<0.05, nominal p<0.10, lowest
       2,000 p-values, lowest 5,000 p-values; split each by direction.
    4. Run ChEA_2022 + TRRUST_Transcription_Factors_2019 enrichment via
       gseapy/Enrichr for each of the 8 pool-direction combinations.
    5. Save gene pools, full per-term enrichment records (with Genes and
       original Term text preserved), summary counts, threshold-stability,
       cross-database agreement, and a run-provenance record.

RUNX2 and SMAD1 still require manual term/overlap auditing before being
placed in the main Wnt-osteogenesis network -- nothing in this script
changes that requirement.

Dependencies:
    pip install gseapy pandas

Usage:
    python 04_threshold_deg_tf_sensitivity.py
"""

from __future__ import annotations

import csv
import datetime
import json
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd

LIBRARIES = ["ChEA_2022", "TRRUST_Transcription_Factors_2019"]
ORGANISM = "mouse"
FDR_CUTOFF = 0.05

CANDIDATE_DESEQ2_PATHS = [
    Path("data/raw/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.txt"),
    Path("data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"),
    Path("data/raw/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"),
]
OUTDIR = Path("results/tf_sensitivity")

SCENARIOS = [
    ("p005", "Nominal p < 0.05", "p", 0.05),
    ("p010", "Nominal p < 0.10", "p", 0.10),
    ("top2000", "Top 2,000 by nominal p", "top", 2000),
    ("top5000", "Top 5,000 by nominal p", "top", 5000),
]

EXPECTED_COUNTS = {
    "p005": {"downregulated": 994, "upregulated": 340},
    "p010": {"downregulated": 1608, "upregulated": 689},
    "top2000": {"downregulated": 1428, "upregulated": 572},
    "top5000": {"downregulated": 3250, "upregulated": 1750},
}


def resolve_deseq2_path() -> Path:
    for p in CANDIDATE_DESEQ2_PATHS:
        if p.exists():
            return p
    raise FileNotFoundError(f"DESeq2 file not found at any of: {[str(p) for p in CANDIDATE_DESEQ2_PATHS]}")


def sniff_delimiter(path: Path) -> str:
    with open(path, "r", newline="") as f:
        sample = f.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        return "\t" if path.suffix == ".txt" else ","


def bh_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(ranked, 0, 1)[np.argsort(order)]


def load_validated_table() -> pd.DataFrame:
    path = resolve_deseq2_path()
    delim = sniff_delimiter(path)
    print(f"Loading {path} (detected delimiter: {'TAB' if delim == chr(9) else repr(delim)})")
    df = pd.read_csv(path, sep=delim)
    df = df.dropna(subset=["pvalue", "symbol", "log2FoldChange"]).copy()
    df["symbol"] = df["symbol"].astype(str).str.strip().str.upper()
    df = df[df.symbol.ne("")]

    df = (
        df.assign(abs_log2fc=df["log2FoldChange"].abs())
          .sort_values(["symbol", "pvalue", "abs_log2fc", "baseMean"], ascending=[True, True, False, False])
          .drop_duplicates("symbol", keep="first")
    )

    # NOTE on scope: BH correction is recomputed AFTER symbol-level
    # deduplication, over this project's established symbol-level audit
    # set -- not over DESeq2's original tested features prior to any
    # deduplication. This is why the column is named padj_recomputed and
    # must never be described as reproducing DESeq2's native padj, which
    # is computed by the algorithm itself over a different (pre-dedup)
    # feature set.
    df["padj_recomputed"] = bh_adjust(df["pvalue"].values)

    n_fdr_sig = int((df["padj_recomputed"] < FDR_CUTOFF).sum())
    print(f"Loaded {len(df)} unique-symbol genes with valid pvalue/log2FoldChange (deterministic dedup).")

    if n_fdr_sig != 30:
        raise RuntimeError(
            f"BH audit produced {n_fdr_sig} significant symbols; expected 30. "
            "Confirm the input file, delimiter, deduplication, and p-value columns "
            "before running enrichment."
        )
    print(f"Recomputed-BH-FDR audit: {n_fdr_sig} genes at padj_recomputed<0.05 -- MATCH (expected 30).")
    return df


def build_gene_pools(df: pd.DataFrame) -> dict:
    pools = {}
    input_count_rows = []
    pool_dir = OUTDIR / "gene_pools"
    pool_dir.mkdir(parents=True, exist_ok=True)

    for scenario, label, mode, value in SCENARIOS:
        if mode == "p":
            selected = df[df.pvalue < value]
        elif mode == "top":
            selected = df.sort_values("pvalue").head(int(value))
        else:
            raise ValueError(f"Unknown selection mode: {mode}")

        down = sorted(selected.loc[selected.log2FoldChange < 0, "symbol"].unique())
        up = sorted(selected.loc[selected.log2FoldChange > 0, "symbol"].unique())
        pools[scenario] = {"downregulated": down, "upregulated": up}

        observed = {"downregulated": len(down), "upregulated": len(up)}
        if observed != EXPECTED_COUNTS[scenario]:
            raise RuntimeError(
                f"{scenario} gene counts differ from the validated workflow: "
                f"observed={observed}, expected={EXPECTED_COUNTS[scenario]}. "
                "Confirm the input file and selection logic before enrichment."
            )
        print(f"  {scenario}: {len(down) + len(up)} total "
              f"(down={len(down)}, up={len(up)}) -- VALIDATED")

        for direction, genes in [("downregulated", down), ("upregulated", up)]:
            input_count_rows.append({
                "scenario": scenario, "scenario_label": label, "selection_mode": mode,
                "selection_value": value, "direction": direction, "gene_count": len(genes),
            })
            pd.DataFrame({"symbol": genes}).to_csv(pool_dir / f"{scenario}_{direction}.csv", index=False)

    counts_df = pd.DataFrame(input_count_rows)
    print("\nInput counts:")
    print(counts_df.to_string(index=False))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    counts_df.to_csv(OUTDIR / "input_counts.csv", index=False)
    return pools


def run_enrichment_for_pool(genes: list, scenario: str, direction: str) -> pd.DataFrame:
    if not genes:
        print(f"  [WARNING] {scenario}/{direction}: empty gene list, skipping.")
        return pd.DataFrame()

    result = gp.enrichr(
        gene_list=genes, gene_sets=LIBRARIES, organism=ORGANISM,
        outdir=None, cutoff=1.0,
    )
    table = result.results.copy()
    table["scenario"] = scenario
    table["direction"] = direction
    table["input_gene_count"] = len(genes)
    return table


def main():
    run_started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    df = load_validated_table()
    pools = build_gene_pools(df)

    all_records = []
    for scenario, directions in pools.items():
        for direction, genes in directions.items():
            print(f"\nRunning enrichment: scenario={scenario}, direction={direction}, n_genes={len(genes)}")
            table = run_enrichment_for_pool(genes, scenario, direction)
            if not table.empty:
                all_records.append(table)

    OUTDIR.mkdir(parents=True, exist_ok=True)

    if not all_records:
        provenance = {
            "run_started": run_started, "libraries": LIBRARIES, "organism": ORGANISM,
            "fdr_cutoff": FDR_CUTOFF, "gseapy_version": gp.__version__,
            "status": "FAILED: no enrichment records returned",
        }
        (OUTDIR / "run_provenance.json").write_text(json.dumps(provenance, indent=2))
        raise RuntimeError(
            "No enrichment records were returned. Check internet connectivity, "
            "Enrichr availability, and library names."
        )

    full = pd.concat(all_records, ignore_index=True)

    full["term_original"] = full["Term"]
    full["regulator"] = full["Term"].str.extract(r"^\s*([A-Za-z0-9_.-]+)", expand=False).str.upper()

    full.to_csv(OUTDIR / "enrichment_full_by_scenario.csv", index=False)
    print(f"\nSaved full per-term enrichment records (with Genes and term_original columns) to "
          f"{OUTDIR / 'enrichment_full_by_scenario.csv'} -- required by tf_target_deg_spider_maps.py.")

    summary_rows = []
    for (scenario, direction, gene_set), group in full.groupby(["scenario", "direction", "Gene_set"]):
        summary_rows.append({
            "scenario": scenario, "direction": direction, "Gene_set": gene_set,
            "tested_records": len(group),
            "significant_records": int((group["Adjusted P-value"] < FDR_CUTOFF).sum()),
            "minimum_adjusted_p": group["Adjusted P-value"].min(),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTDIR / "enrichment_summary.csv", index=False)
    print("\nEnrichment summary:")
    print(summary_df.to_string(index=False))

    stability_rows = []
    n_scenarios_total = len(SCENARIOS)
    for (regulator, direction, gene_set), group in full.groupby(["regulator", "direction", "Gene_set"]):
        sig_group = group[group["Adjusted P-value"] < FDR_CUTOFF]
        scenarios_sig = sorted(sig_group.scenario.unique())
        stability_rows.append({
            "direction": direction, "Gene_set": gene_set, "regulator": regulator,
            "scenarios_tested": group.scenario.nunique(),
            "scenarios_FDR_significant": len(scenarios_sig),
            "minimum_adjusted_p": group["Adjusted P-value"].min(),
            "significant_scenarios": ";".join(scenarios_sig) if scenarios_sig else None,
            "stable_all_four": len(scenarios_sig) == n_scenarios_total,
        })
    stability_df = pd.DataFrame(stability_rows).sort_values(
        ["scenarios_FDR_significant", "minimum_adjusted_p"], ascending=[False, True])
    stability_df.to_csv(OUTDIR / "regulator_threshold_stability.csv", index=False)
    print(f"\nSaved regulator threshold stability ({len(stability_df)} regulator/direction/library rows)")

    agreement_columns = ["scenario", "direction", "regulator", "significant_libraries", "significant_in_both_databases"]
    agreement_rows = []
    for (scenario, direction, regulator), group in full.groupby(["scenario", "direction", "regulator"]):
        sig_libs = sorted(group.loc[group["Adjusted P-value"] < FDR_CUTOFF, "Gene_set"].unique())
        if len(sig_libs) >= 2:
            agreement_rows.append({
                "scenario": scenario, "direction": direction, "regulator": regulator,
                "significant_libraries": ";".join(sig_libs), "significant_in_both_databases": True,
            })
    agreement_df = pd.DataFrame(agreement_rows, columns=agreement_columns)
    agreement_df.to_csv(OUTDIR / "cross_database_agreement_by_scenario.csv", index=False)
    print(f"\nDual-database scenario-direction-regulator hits: {len(agreement_df)}")
    if not agreement_df.empty:
        print(agreement_df.to_string(index=False))
    print("NOTE: RUNX2 and SMAD1 (if present above) still require manual term/overlap "
          "auditing before being placed in the main Wnt-osteogenesis network.")

    provenance = {
        "run_started": run_started,
        "run_completed": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "libraries": LIBRARIES,
        "organism": ORGANISM,
        "fdr_cutoff": FDR_CUTOFF,
        "gseapy_version": gp.__version__,
        "input_gene_counts_by_scenario_direction": {
            f"{s}_{d}": len(pools[s][d]) for s in pools for d in pools[s]
        },
        "total_enrichment_records_returned": len(full),
        "status": "OK",
    }
    (OUTDIR / "run_provenance.json").write_text(json.dumps(provenance, indent=2))
    print(f"\nSaved run provenance to {OUTDIR / 'run_provenance.json'}")
    print(f"All outputs saved to {OUTDIR}")


if __name__ == "__main__":
    main()
