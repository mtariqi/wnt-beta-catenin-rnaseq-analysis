import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_canonical_input_schema_and_counts():
    path = ROOT / "data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
    table = pd.read_csv(path)
    assert list(table.columns) == [
        "gene_id",
        "symbol",
        "baseMean",
        "log2FoldChange",
        "pvalue",
    ]
    assert len(table) == 22_444
    valid = table.dropna(subset=["symbol", "log2FoldChange", "pvalue"])
    assert len(valid) == 21_730
    assert valid["symbol"].nunique() == 21_715


def test_deg_outputs():
    result_dir = ROOT / "results/deg"
    all_degs = pd.read_csv(result_dir / "all_DEGs.csv")
    up = pd.read_csv(result_dir / "upregulated_genes.csv")
    down = pd.read_csv(result_dir / "downregulated_genes.csv")
    strict = pd.read_csv(result_dir / "strict_DEGs_p0.01_abslog2fc_gt2.csv")

    assert (len(all_degs), len(up), len(down)) == (745, 84, 661)
    assert len(up) + len(down) == len(all_degs)
    assert len(strict) == 249
    assert (strict["pvalue"] < 0.01).all()
    assert (strict["log2FoldChange"].abs() > 2).all()
    assert not strict["symbol"].duplicated().any()


def test_custom_wnt_signature_outputs():
    summary = pd.read_csv(
        ROOT / "results/wnt_signatures/wnt_signature_summary.csv"
    )
    assert len(summary) == 4
    assert summary["fisher_FDR_BH"].between(0, 1).all()
    assert (
        summary["upregulated_overlap"] + summary["downregulated_overlap"]
        == summary["deg_overlap"]
    ).all()

    beta = summary.loc[summary["signature"] == "beta_catenin_dependent"].iloc[0]
    tcf = summary.loc[
        summary["signature"].isin(["tcf_dependent", "tcf_independent"])
    ]
    assert beta["detectable_in_dataset"] == tcf["detectable_in_dataset"].sum()
    assert beta["deg_overlap"] == tcf["deg_overlap"].sum()


def test_custom_preranked_gsea():
    result_dir = ROOT / "results/gsea"
    ranked = pd.read_csv(
        result_dir / "full_transcriptome_ranked.rnk",
        sep="\t",
        header=None,
        names=["symbol", "score"],
        float_precision="round_trip",
    )
    results = pd.read_csv(result_dir / "custom_wnt_preranked_gsea.csv")

    assert len(ranked) == 21_715
    assert not ranked["symbol"].duplicated().any()
    assert not ranked["score"].duplicated().any()
    assert ranked["score"].is_monotonic_decreasing

    epithelial = results.loc[
        results["Term"] == "EPITHELIAL_JUNCTION_CORE"
    ].iloc[0]
    assert epithelial["NES"] < 0
    assert epithelial["FDR q-val"] < 0.05

    literature_terms = [
        "BetaCatenin_Independent",
        "BetaCatenin_Dependent",
        "TCF_Independent",
        "TCF_Dependent",
    ]
    literature = results.loc[results["Term"].isin(literature_terms)]
    assert len(literature) == 4
    assert (literature["FDR q-val"] >= 0.05).all()


def test_broad_pathway_gsea_negative_result():
    result_dir = ROOT / "results/pathway_gsea"
    combined = pd.read_csv(result_dir / "combined_pathway_preranked_gsea.csv")
    significant = pd.read_csv(result_dir / "significant_pathways_FDR005.csv")

    assert len(combined) == 5_057
    assert significant.empty
    assert set(combined["Library"]) == {
        "GO_BP_2023",
        "Reactome_2022",
        "KEGG_2026",
    }
    assert set(combined["FDR_scope"]) == {
        "within_GO_BP_2023",
        "within_Reactome_2022",
        "within_KEGG_2026",
    }


def test_craniofacial_outputs():
    result_dir = ROOT / "results/craniofacial_skeletal"
    terms = pd.read_csv(result_dir / "01_selected_craniofacial_skeletal_GO_terms.csv")
    universe = pd.read_csv(
        result_dir / "03_systematic_craniofacial_skeletal_gene_universe.csv"
    )
    intersection = pd.read_csv(
        result_dir / "04_craniofacial_skeletal_DESeq2_intersection.csv"
    )
    references = pd.read_csv(result_dir / "07_canonical_Wnt_reference_genes.csv")

    assert terms["GO_term"].nunique() == 44
    assert len(universe) == 496
    assert len(intersection) == 453
    assert intersection["symbol"].nunique() == 453

    fdr_genes = set(
        intersection.loc[intersection["BH_FDR_AUDIT"] < 0.05, "symbol"]
    )
    assert fdr_genes == {"Bmp6", "Lama5"}
    assert len(references) == 6
    assert (references["Significance"] == "Not significant").all()


def test_upstream_tf_target_enrichment():
    result_dir = ROOT / "results/enrichment"
    metadata = json.loads((result_dir / "enrichr_run_metadata.json").read_text())
    down = pd.read_csv(result_dir / "upstream_downregulated_FDR05.csv")
    up = pd.read_csv(result_dir / "upstream_upregulated_FDR05.csv")
    plotted = pd.read_csv(result_dir / "verified_slide10_plotted_values.csv")
    agreement = pd.read_csv(result_dir / "cross_database_agreement.csv")

    assert metadata["input_definition"] == "nominal p < 0.05; no fold-change cutoff"
    assert metadata["downregulated_input_genes"] == 994
    assert metadata["upregulated_input_genes"] == 340
    assert len(down) == 72
    assert len(up) == 5
    assert (plotted["Adjusted P-value"] < 0.05).all()
    assert agreement["significant_in_both"].sum() == 0
