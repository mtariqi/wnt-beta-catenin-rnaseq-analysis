#!/usr/bin/env python3
"""Append audited TF sensitivity results to a report and write a Professor email."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report", required=True, type=Path, help="Existing validated .docx report")
    p.add_argument("--results", type=Path, default=Path("results/tf_sensitivity"))
    p.add_argument("--output-report", type=Path, default=Path("reports/manuscripts/Wnt_TF_sensitivity_report_updated.docx"))
    p.add_argument("--email", type=Path, default=Path("reports/Professor_Hsu_TF_sensitivity_email.txt"))
    return p.parse_args()


def add_df_table(doc: Document, df: pd.DataFrame, columns: list[str]) -> None:
    table = doc.add_table(rows=1, cols=len(columns))
    table.style = "Table Grid"
    for i, col in enumerate(columns):
        table.rows[0].cells[i].text = col
    for _, row in df.iterrows():
        cells = table.add_row().cells
        for i, col in enumerate(columns):
            value = row[col]
            if isinstance(value, float):
                cells[i].text = f"{value:.3g}"
            else:
                cells[i].text = str(value)
    for row in table.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(8)


def main() -> None:
    a = parse_args()
    if not a.report.exists():
        raise FileNotFoundError(a.report)
    counts = pd.read_csv(a.results / "threshold_input_counts.csv")
    dual = pd.read_csv(a.results / "dual_database_regulators_by_scenario.csv")
    details = pd.read_csv(a.results / "dual_database_significant_details.csv")
    stability = pd.read_csv(a.results / "audited_regulator_threshold_stability.csv")
    heatmap = a.results / "tf_threshold_stability_heatmap.png"
    if not heatmap.exists():
        raise FileNotFoundError("Run script 15 before building the report.")

    total_by_scenario = counts.groupby(["scenario", "scenario_label"], as_index=False)["gene_count"].sum()
    stable = stability[stability["stable_all_four"]].copy()
    dual_names = sorted(dual["regulator"].unique()) if len(dual) else []

    doc = Document(a.report)
    doc.add_page_break()
    doc.add_heading("Supplementary analysis: TF enrichment threshold sensitivity", level=1)
    doc.add_paragraph(
        "To evaluate dependence on the input-gene definition, upstream transcription-factor target-set enrichment "
        "was repeated using nominal p < 0.05, nominal p < 0.10, the 2,000 smallest nominal p-values and the 5,000 "
        "smallest nominal p-values. No absolute log2-fold-change cutoff was applied. Upregulated and downregulated "
        "genes were submitted separately to ChEA 2022 and TRRUST Transcription Factors 2019."
    )

    doc.add_heading("Input-list sensitivity", level=2)
    add_df_table(doc, total_by_scenario, ["scenario_label", "gene_count"])
    doc.add_paragraph(
        "The nominal p < 0.10 definition included 2,297 genes and therefore directly addressed the request for an "
        "input exceeding 2,000 genes. The top-5,000 analysis evaluated a substantially broader input universe."
    )

    doc.add_heading("Cross-threshold and cross-database findings", level=2)
    doc.add_paragraph(
        f"The audited analysis identified {len(dual)} scenarioâ€“directionâ€“regulator combinations that met term-level "
        f"BH-adjusted p < 0.05 in both databases. The implicated parsed regulator labels were: "
        f"{', '.join(dual_names) if dual_names else 'none'}. Cross-database agreement was threshold-dependent and "
        "must not be generalized to the original nominal p < 0.05 analysis unless present in that scenario."
    )
    if len(dual):
        show = dual[["scenario", "direction", "regulator", "significant_databases"]].copy()
        add_df_table(doc, show, list(show.columns))

    doc.add_heading("Three independent evidence layers", level=2)
    layers = pd.DataFrame([
        ["Cross-threshold stability", "Target-set enrichment recurs across input definitions", "Does not prove TF activity"],
        ["Cross-database agreement", "ChEA and TRRUST both support the association", "Databases are not biologically independent experiments"],
        ["Literature verification", "Published evidence supports a proposed pathway relationship", "Independent of this dataset's enrichment statistics"],
    ], columns=["Evidence layer", "What it supports", "What it does not establish"])
    add_df_table(doc, layers, list(layers.columns))

    doc.add_heading("RUNX2 and SMAD1 interpretation", level=2)
    rs = details[details["regulator"].isin(["RUNX2", "SMAD1"])].copy()
    doc.add_paragraph(
        "RUNX2- and SMAD1-associated target sets are retained as enrichment findings only. They were not added to "
        "the literature network because the exact database terms, overlap genes, enrichment statistics and "
        "citation-supported pathway edges require a separate biological audit. Target-set enrichment does not "
        "demonstrate increased TF expression, phosphorylation, nuclear localization, direct binding or activity."
    )
    if len(rs):
        cols = [c for c in ["scenario", "direction", "regulator", "Gene_set", "Term", "Overlap", "Adjusted P-value"] if c in rs]
        add_df_table(doc, rs[cols], cols)

    doc.add_heading("Threshold-stability visualization", level=2)
    paragraph = doc.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.add_run().add_picture(str(heatmap), width=Inches(6.5))
    caption = doc.add_paragraph(
        "Figure S1. Threshold stability of transcription-factor target-set enrichment. Color represents "
        "âˆ’log10 of the BH-adjusted enrichment p-value; filled circles identify FDR < 0.05. Columns distinguish "
        "the input definition, gene-list direction and annotation database. The heatmap compares sensitivity "
        "but does not constitute an independent inferential analysis."
    )
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.add_heading("Interpretive conclusion", level=2)
    doc.add_paragraph(
        f"Several TF target-set associations were stable across all four input definitions ({len(stable)} "
        "directionâ€“databaseâ€“regulator combinations). Broader input definitions revealed additional cross-database "
        "agreement, but the number and identity of significant terms varied by database and threshold. Consequently, "
        "robust associations should be prioritized for experimental follow-up while retaining explicit database, "
        "direction and threshold qualifiers."
    )
    a.output_report.parent.mkdir(parents=True, exist_ok=True)
    doc.save(a.output_report)

    p005_dual = dual[dual["scenario"].eq("p005")]
    broader_dual = dual[~dual["scenario"].eq("p005")]
    email = f"""Subject: TF target-set enrichment sensitivity and cross-database robustness update

Dear Professor Hsu,

I completed the requested sensitivity analysis using four input definitions: nominal p < 0.05, nominal p < 0.10, the 2,000 smallest nominal p-values, and the 5,000 smallest nominal p-values, without an absolute fold-change cutoff. These definitions contained 1,334, 2,297, 2,000, and 5,000 unique gene symbols, respectively.

The original nominal p < 0.05 result remains correctly scoped: it contained {len(p005_dual)} scenario-direction-regulator combinations significant in both ChEA 2022 and TRRUST 2019. The broader input definitions yielded {len(broader_dual)} additional cross-database combinations. Across all scenarios, {len(dual)} dual-database combinations were identified, involving: {', '.join(dual_names) if dual_names else 'none'}.

Several associations also remained significant across all four input definitions, including stable Wnt-related and osteogenic target-set patterns. I have separated three forms of support in the updated report: cross-threshold stability, cross-database agreement, and literature verification. These layers are complementary but not interchangeable.

RUNX2- and SMAD1-associated target sets are reported as enrichment findings only. I have not added them to the literature network pending a dedicated audit of their exact database terms, overlapping genes, enrichment FDR values, and citation-supported pathway relationships. Likewise, enrichment does not establish TF expression, activation, direct binding, regulatory direction, or causality.

I have attached the updated report, audited result tables, and threshold-stability heatmap for review.

Best regards,
Tariqul Islam
"""
    a.email.parent.mkdir(parents=True, exist_ok=True)
    a.email.write_text(email, encoding="utf-8")
    print(f"Saved report: {a.output_report}")
    print(f"Saved email: {a.email}")


if __name__ == "__main__":
    main()
