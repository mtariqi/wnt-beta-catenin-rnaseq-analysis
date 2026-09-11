import numpy as np
import pandas as pd
import gseapy as gp

from pathlib import Path

# ============================================================
# Paths
# ============================================================

DESEQ_FILE = Path(
    "data/processed/"
    "DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
)

OUTPUT_DIR = Path("results/pathway_gsea")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load DESeq2 summary results
# ============================================================

df = pd.read_csv(DESEQ_FILE)
df.columns = df.columns.str.strip()

df = df.dropna(
    subset=["symbol", "log2FoldChange", "pvalue"]
).copy()

df["symbol_upper"] = (
    df["symbol"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# ============================================================
# Exploratory ranking statistic
#
# log2FC * -log10(p)
#
# Positive = mutant-upregulated side
# Negative = mutant-downregulated side
# ============================================================

tiny = np.finfo(float).tiny

df["rank_score"] = (
    df["log2FoldChange"]
    * -np.log10(df["pvalue"].clip(lower=tiny))
)

# Resolve duplicate gene symbols by retaining the strongest
# absolute ranking statistic
ranked = (
    df.assign(abs_rank=df["rank_score"].abs())
    .sort_values("abs_rank", ascending=False)
    .drop_duplicates("symbol_upper")
    [["symbol_upper", "rank_score"]]
    .sort_values("rank_score", ascending=False)
)

# Deterministically resolve exact ties caused by rounded source statistics.
ranked = ranked.sort_values(["rank_score", "symbol_upper"], ascending=[False, True])
tie_index = ranked.groupby("rank_score", sort=False).cumcount().to_numpy()
ranked["rank_score"] = (
    ranked["rank_score"].to_numpy()
    - tie_index * np.spacing(np.abs(ranked["rank_score"].to_numpy()))
)
ranked = ranked.sort_values("rank_score", ascending=False)

print(f"Genes in ranked transcriptome: {len(ranked):,}")

ranked_file = OUTPUT_DIR / "full_transcriptome_ranked.rnk"

ranked.to_csv(
    ranked_file,
    sep="\t",
    index=False,
    header=False
)

print(f"Saved ranking: {ranked_file}")

# ============================================================
# Mouse pathway libraries
# ============================================================

libraries = {
    "GO_BP_2023": "GO_Biological_Process_2023",
    "Reactome_2022": "Reactome_2022",
    "KEGG_2026": "KEGG_2026",
}

all_results = []

# ============================================================
# Run GSEA library by library
# ============================================================

for short_name, library_name in libraries.items():

    print("\n" + "=" * 70)
    print(f"Loading {library_name}")
    print("=" * 70)

    gene_sets = gp.get_library(
        name=library_name,
        organism="mouse"
    )

    print(
        f"{short_name}: {len(gene_sets):,} gene sets loaded"
    )

    library_snapshot = pd.DataFrame(
        [(term, gene) for term, genes in gene_sets.items() for gene in genes],
        columns=["Term", "Gene"],
    )
    library_snapshot.to_csv(
        OUTPUT_DIR / f"{short_name}_gene_set_membership.csv",
        index=False,
    )

    pre_res = gp.prerank(
        rnk=ranked,
        gene_sets=gene_sets,
        min_size=10,
        max_size=500,
        permutation_num=1000,
        weight=1.0,
        threads=1,
        seed=42,
        outdir=None,
        verbose=False
    )

    results = pre_res.res2d.copy()

    results["Library"] = short_name
    results["FDR_scope"] = f"within_{short_name}"

    output_file = (
        OUTPUT_DIR /
        f"{short_name}_preranked_gsea.csv"
    )

    results.to_csv(
        output_file,
        index=False
    )

    print(f"Saved: {output_file}")

    all_results.append(results)

# ============================================================
# Combine results
# ============================================================

combined = pd.concat(
    all_results,
    ignore_index=True
)

combined_file = (
    OUTPUT_DIR /
    "combined_pathway_preranked_gsea.csv"
)

combined.to_csv(
    combined_file,
    index=False
)

print(f"\nSaved combined results: {combined_file}")

# ============================================================
# Significant pathways
# ============================================================

sig = combined[
    combined["FDR q-val"] < 0.05
].copy()

sig = sig.sort_values(
    ["FDR q-val", "NES"]
)

sig_file = (
    OUTPUT_DIR /
    "significant_pathways_FDR005.csv"
)

sig.to_csv(
    sig_file,
    index=False
)

print(
    f"Significant pathways (FDR < 0.05): {len(sig):,}"
)

# ============================================================
# Targeted biological themes
# ============================================================

keywords = [
    "WNT",
    "BETA CATENIN",
    "BETA-CATENIN",
    "EPITHEL",
    "JUNCTION",
    "ADHESION",
    "DESMOSOME",
    "KERATIN",
    "CYTOSKELETON",
    "POLARITY",
    "PLANAR CELL",
    "INTERFERON",
    "CYTOKINE",
    "INFLAM",
    "FOXO",
    "OXIDATIVE",
    "STRESS",
    "DIFFERENTIATION",
]

pattern = "|".join(keywords)

targeted = combined[
    combined["Term"]
    .astype(str)
    .str.upper()
    .str.contains(pattern, regex=True, na=False)
].copy()

targeted = targeted.sort_values(
    ["FDR q-val", "NES"]
)

targeted_file = (
    OUTPUT_DIR /
    "targeted_biological_pathways.csv"
)

targeted.to_csv(
    targeted_file,
    index=False
)

print(
    f"Targeted biological pathways: {len(targeted):,}"
)

# ============================================================
# Print strongest results
# ============================================================

show_cols = [
    c for c in [
        "Library",
        "Term",
        "NES",
        "NOM p-val",
        "FDR q-val",
        "Lead_genes"
    ]
    if c in combined.columns
]

print("\n===== TOP FDR-SIGNIFICANT PATHWAYS =====\n")

if len(sig) > 0:

    print(
        sig[
            show_cols
        ].head(30).to_string(index=False)
    )

else:
    print("No pathways reached FDR < 0.05.")

print("\n===== TARGETED BIOLOGICAL THEMES =====\n")

print(
    targeted[
        show_cols
    ].head(40).to_string(index=False)
)
