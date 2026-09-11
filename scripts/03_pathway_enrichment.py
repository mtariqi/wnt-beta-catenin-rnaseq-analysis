import pandas as pd
import gseapy as gp
from pathlib import Path

UP_FILE = Path("results/deg/upregulated_gene_symbols.txt")
DOWN_FILE = Path("results/deg/downregulated_gene_symbols.txt")

OUTPUT_DIR = Path("results/pathway_enrichment")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_genes(path):
    genes = (
        pd.read_csv(path, header=None)[0]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    return list(dict.fromkeys(genes))

up_genes = load_genes(UP_FILE)
down_genes = load_genes(DOWN_FILE)

pd.Series(up_genes, name="symbol").to_csv(OUTPUT_DIR / "upregulated_input_genes.csv", index=False)
pd.Series(down_genes, name="symbol").to_csv(OUTPUT_DIR / "downregulated_input_genes.csv", index=False)

print(f"Upregulated genes: {len(up_genes)}")
print(f"Downregulated genes: {len(down_genes)}")

gene_sets = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "KEGG_2026",
    "Reactome_2022",
]

def run_enrichment(genes, label):
    print(f"\n=== {label.upper()} ENRICHMENT ===")

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism="mouse",
        outdir=None,
        cutoff=1.0,
        no_plot=True
    )

    results = enr.results.copy()

    output_file = OUTPUT_DIR / f"{label}_all_results.csv"
    results.to_csv(output_file, index=False)

    print(f"Saved: {output_file}")

    if results.empty:
        print("No enrichment terms returned.")
        return

    # Sort by adjusted p-value if available
    if "Adjusted P-value" in results.columns:
        results = results.sort_values("Adjusted P-value")

        significant = results[results["Adjusted P-value"] < 0.05].copy()
        significant.to_csv(OUTPUT_DIR / f"{label}_FDR05.csv", index=False)
        print(f"FDR < 0.05 records: {len(significant)}")

    show_cols = [
        c for c in [
            "Gene_set",
            "Term",
            "Adjusted P-value",
            "P-value",
            "Odds Ratio",
            "Combined Score",
            "Genes"
        ]
        if c in results.columns
    ]

    print(results[show_cols].head(20).to_string(index=False))


run_enrichment(up_genes, "upregulated")
run_enrichment(down_genes, "downregulated")
