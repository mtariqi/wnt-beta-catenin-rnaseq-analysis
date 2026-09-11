#!/usr/bin/env python3
"""Reproducible TRRUST/ChEA target-set enrichment for the Wnt mouse DEG study.

The script:
  1. reads a DESeq2 results table;
  2. selects nominal p < 0.05 genes without a fold-change magnitude cutoff;
  3. validates the requested Enrichr libraries;
  4. runs separate up/down enrichment analyses;
  5. exports complete and FDR < 0.05 results;
  6. audits cross-database agreement;
  7. plots verified mouse-labelled records; and
  8. optionally replaces the final slide of a PowerPoint deck.

Important: gp.enrichr() submits gene symbols to the Enrichr web service.

Dependencies:
    pip install pandas numpy matplotlib gseapy python-pptx

Example:
    python 13_upstream_tf_target_enrichment.py \
      --input DESeq2_results.csv \
      --outdir results/enrichment \
      --pptx Sep09_presentation_revised_professor_feedback.pptx \
      --slide-number 10 \
      --output-pptx Sep09_presentation_verified_slide10.pptx
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


LIBRARIES = ["TRRUST_Transcription_Factors_2019", "ChEA_2022"]
LIBRARY_COLORS = {
    "ChEA_2022": "#2176C7",
    "TRRUST_Transcription_Factors_2019": "#E66B20",
}
REQUIRED_COLUMNS = {"symbol", "pvalue", "log2FoldChange"}
ENRICHR_COLUMNS = {
    "Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Odds Ratio"
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run and audit TF-target enrichment for up/down mouse DEGs."
    )
    parser.add_argument("--input", required=True, type=Path, help="DESeq2 CSV/TSV/TXT file")
    parser.add_argument("--outdir", type=Path, default=Path("results/enrichment"))
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--top-down", type=int, default=8)
    parser.add_argument("--organism", default="mouse")
    parser.add_argument("--pptx", type=Path, help="Optional deck containing the slide to replace")
    parser.add_argument(
        "--slide-number",
        type=int,
        help="1-based slide number to replace; required when --pptx is supplied",
    )
    parser.add_argument("--output-pptx", type=Path, help="Required when --pptx is supplied")
    return parser.parse_args()


def read_deseq2(path: Path) -> tuple[pd.DataFrame, int]:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    table = pd.read_csv(path, sep=sep)
    missing = REQUIRED_COLUMNS - set(table.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    table = table.dropna(subset=list(REQUIRED_COLUMNS)).copy()
    valid_rows = len(table)
    table["symbol"] = table["symbol"].astype(str).str.strip().str.upper()
    table = table[table.symbol.ne("")]
    # Make duplicate resolution explicit: retain the row with the smallest
    # nominal p-value for each gene symbol, independent of source-file order.
    table = table.sort_values("pvalue", kind="stable").drop_duplicates("symbol")
    return table, valid_rows


def select_inputs(table: pd.DataFrame, p_cutoff: float) -> tuple[list[str], list[str]]:
    significant = table[table.pvalue < p_cutoff]
    up = sorted(significant.loc[significant.log2FoldChange > 0, "symbol"].unique())
    down = sorted(significant.loc[significant.log2FoldChange < 0, "symbol"].unique())
    if not up or not down:
        raise ValueError("The selected threshold produced an empty up- or downregulated list.")
    return up, down


def validate_libraries(organism: str) -> dict[str, bool]:
    available = gp.get_library_name(organism=organism)
    status = {name: name in available for name in LIBRARIES}
    if not all(status.values()):
        raise RuntimeError(f"Enrichr library validation failed: {status}")
    return status


def library_coverage(
    up: list[str], down: list[str], organism: str
) -> pd.DataFrame:
    records = []
    for library in LIBRARIES:
        definitions = gp.get_library(name=library, organism=organism)
        library_symbols = {
            str(gene).strip().upper()
            for members in definitions.values()
            for gene in members
            if str(gene).strip()
        }
        for direction, genes in (("upregulated", up), ("downregulated", down)):
            overlap = set(genes) & library_symbols
            records.append({
                "library": library,
                "direction": direction,
                "input_genes": len(genes),
                "input_genes_present_in_library": len(overlap),
                "coverage_percent": 100 * len(overlap) / len(genes),
            })
    return pd.DataFrame(records)


def run_enrichr(
    genes: list[str], direction: str, organism: str, outdir: Path, fdr_cutoff: float
) -> tuple[pd.DataFrame, pd.DataFrame]:
    result = gp.enrichr(
        gene_list=genes,
        gene_sets=LIBRARIES,
        organism=organism,
        outdir=str(outdir / f"upstream_{direction}"),
        cutoff=1.0,
        no_plot=True,
    )
    table = result.results.copy()
    if table.empty:
        raise RuntimeError(f"Enrichr returned no results for {direction} genes.")
    missing = ENRICHR_COLUMNS - set(table.columns)
    if missing:
        raise RuntimeError(f"Enrichr output is missing columns: {sorted(missing)}")

    table["DEG_direction"] = direction
    table["Input_gene_count"] = len(genes)
    table = table.sort_values(["Adjusted P-value", "P-value", "Gene_set", "Term"])
    significant = table[table["Adjusted P-value"] < fdr_cutoff].copy()

    table.to_csv(outdir / f"upstream_{direction}_all_results.csv", index=False)
    significant.to_csv(outdir / f"upstream_{direction}_FDR05.csv", index=False)
    return table, significant


def regulator_from_term(term: str) -> str:
    match = re.match(r"^([^ _]+)", str(term))
    if not match:
        raise ValueError(f"Cannot parse regulator from term: {term}")
    return match.group(1)


def build_agreement_table(
    down_fdr: pd.DataFrame, up_fdr: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for direction, table in (("downregulated", down_fdr), ("upregulated", up_fdr)):
        work = table.copy()
        work["regulator"] = work.Term.map(regulator_from_term)
        for regulator, group in work.groupby("regulator"):
            libraries = sorted(group.Gene_set.unique())
            rows.append({
                "direction": direction,
                "regulator": regulator,
                "significant_in_both": set(LIBRARIES).issubset(libraries),
                "significant_libraries": "; ".join(libraries),
                "minimum_adjusted_p_value": group["Adjusted P-value"].min(),
            })
    return pd.DataFrame(rows).sort_values(
        ["direction", "significant_in_both", "minimum_adjusted_p_value"],
        ascending=[True, False, True],
    )


def select_plot_records(table: pd.DataFrame, limit: int | None = None) -> pd.DataFrame:
    # The mouse Enrichr endpoint can return mixed-species records. Presentation
    # results are restricted to terms explicitly labelled Mouse.
    selected = table[
        table.Term.astype(str).str.contains(r"\bmouse\b", case=False, regex=True)
    ].copy()
    selected["Regulator"] = selected.Term.map(regulator_from_term)
    selected = (
        selected.sort_values(["Adjusted P-value", "P-value"])
        .drop_duplicates("Regulator")
    )
    selected["bar_value"] = -np.log10(selected["Adjusted P-value"])
    return selected.head(limit) if limit else selected


def plot_results(
    down: pd.DataFrame,
    up: pd.DataFrame,
    output: Path,
    p_cutoff: float,
) -> pd.DataFrame:
    fig, axes = plt.subplots(
        1, 2, figsize=(13.333, 7.5), gridspec_kw={"wspace": 0.38}
    )
    for ax, data, title in (
        (axes[0], down, "Downregulated DEG targets"),
        (axes[1], up, "Upregulated DEG targets"),
    ):
        plotting = data.sort_values("bar_value")
        bars = ax.barh(
            plotting.Regulator,
            plotting.bar_value,
            color=[LIBRARY_COLORS[x] for x in plotting.Gene_set],
            height=0.66,
        )
        ax.axvline(-np.log10(0.05), color="#667085", linestyle="--", linewidth=1.2)
        ax.set_title(title, fontsize=17, fontweight="bold", color="#102A56", pad=12)
        ax.set_xlabel(
            "−log$_{10}$(BH-adjusted enrichment p-value)",
            fontsize=11.5,
            fontweight="bold",
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=11)
        xmax = plotting.bar_value.max() * 1.26
        ax.set_xlim(0, xmax)
        for bar, (_, row) in zip(bars, plotting.iterrows()):
            ax.text(
                bar.get_width() + xmax * 0.018,
                bar.get_y() + bar.get_height() / 2,
                f"FDR={row['Adjusted P-value']:.2g}",
                va="center",
                fontsize=9,
                color="#344054",
            )

    fig.suptitle(
        "Exploratory enrichment of transcription-factor target sets",
        fontsize=23,
        fontweight="bold",
        color="#102A56",
        y=0.965,
    )
    fig.text(
        0.5, 0.905,
        "Mouse-labelled records only; one bar per TF selected by its smallest BH-adjusted p-value",
        ha="center", fontsize=11, color="#475467",
    )
    fig.legend(
        handles=[
            Patch(color=LIBRARY_COLORS["ChEA_2022"], label="ChEA 2022"),
            Patch(
                color=LIBRARY_COLORS["TRRUST_Transcription_Factors_2019"],
                label="TRRUST 2019",
            ),
        ],
        loc="lower center", bbox_to_anchor=(0.5, 0.067), ncol=2,
        frameon=False, fontsize=10.5,
    )
    fig.text(
        0.5, 0.027,
        f"Inputs: {down.Input_gene_count.iloc[0]} downregulated and "
        f"{up.Input_gene_count.iloc[0]} upregulated genes at nominal p < {p_cutoff:g}; "
        "no fold-change cutoff. Target-set enrichment does not establish regulator activation or inhibition.",
        ha="center", fontsize=9.2, color="#667085",
    )
    fig.subplots_adjust(left=0.075, right=0.965, top=0.83, bottom=0.16)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)

    return pd.concat([
        down.assign(Direction="downregulated"),
        up.assign(Direction="upregulated"),
    ])[[
        "Direction", "Regulator", "Gene_set", "Term", "Overlap", "P-value",
        "Adjusted P-value", "Odds Ratio", "bar_value",
    ]]


def replace_slide(
    source: Path, figure: Path, destination: Path, slide_number: int
) -> None:
    try:
        from pptx import Presentation
    except ImportError as exc:
        raise RuntimeError("Install python-pptx to update a presentation.") from exc

    presentation = Presentation(source)
    if not presentation.slides:
        raise ValueError("The supplied PowerPoint contains no slides.")
    if not 1 <= slide_number <= len(presentation.slides):
        raise ValueError(
            f"--slide-number must be between 1 and {len(presentation.slides)}; "
            f"received {slide_number}."
        )
    slide = presentation.slides[slide_number - 1]
    for shape in list(slide.shapes):
        slide.shapes._spTree.remove(shape._element)
    slide.shapes.add_picture(
        str(figure), 0, 0,
        width=presentation.slide_width,
        height=presentation.slide_height,
    )
    presentation.save(destination)


def main() -> None:
    args = parse_args()
    if args.pptx and (not args.output_pptx or args.slide_number is None):
        raise ValueError(
            "--output-pptx and --slide-number are required when --pptx is supplied."
        )
    if not args.pptx and (args.output_pptx or args.slide_number is not None):
        raise ValueError("--pptx is required with --output-pptx or --slide-number.")
    args.outdir.mkdir(parents=True, exist_ok=True)

    table, valid_rows = read_deseq2(args.input)
    up_genes, down_genes = select_inputs(table, args.p_cutoff)
    pd.Series(up_genes, name="symbol").to_csv(
        args.outdir / "upregulated_input_genes.csv", index=False
    )
    pd.Series(down_genes, name="symbol").to_csv(
        args.outdir / "downregulated_input_genes.csv", index=False
    )

    registry = validate_libraries(args.organism)
    coverage = library_coverage(up_genes, down_genes, args.organism)
    coverage.to_csv(args.outdir / "library_symbol_coverage.csv", index=False)

    _, down_fdr = run_enrichr(
        down_genes, "downregulated", args.organism, args.outdir, args.fdr_cutoff
    )
    _, up_fdr = run_enrichr(
        up_genes, "upregulated", args.organism, args.outdir, args.fdr_cutoff
    )

    agreement = build_agreement_table(down_fdr, up_fdr)
    agreement.to_csv(args.outdir / "cross_database_agreement.csv", index=False)

    down_plot = select_plot_records(down_fdr, args.top_down)
    up_plot = select_plot_records(up_fdr)
    figure = args.outdir / "verified_tf_target_enrichment.png"
    plotted = plot_results(down_plot, up_plot, figure, args.p_cutoff)
    plotted.to_csv(args.outdir / "verified_slide10_plotted_values.csv", index=False)

    if args.pptx:
        replace_slide(args.pptx, figure, args.output_pptx, args.slide_number)

    metadata = {
        "python_version": sys.version.split()[0],
        "gseapy_version": gp.__version__,
        "organism": args.organism,
        "libraries": LIBRARIES,
        "library_registry_check": registry,
        "input_definition": f"nominal p < {args.p_cutoff}; no fold-change cutoff",
        "valid_input_rows": valid_rows,
        "unique_analyzed_gene_symbols": len(table),
        "upregulated_input_genes": len(up_genes),
        "downregulated_input_genes": len(down_genes),
        "downregulated_FDR05_records": len(down_fdr),
        "upregulated_FDR05_records": len(up_fdr),
        "regulators_significant_in_both_databases": int(agreement.significant_in_both.sum()),
        "presentation_records_restricted_to_mouse_label": True,
    }
    (args.outdir / "enrichr_run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    (args.outdir / "gseapy_version.txt").write_text(
        gp.__version__ + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    print(f"\nFigure: {figure}")
    if args.output_pptx:
        print(f"Presentation: {args.output_pptx}")


if __name__ == "__main__":
    main()
