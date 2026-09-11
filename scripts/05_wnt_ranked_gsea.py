import numpy as np
import pandas as pd
import gseapy as gp

from pathlib import Path

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

DESEQ_FILE = Path(
    "data/processed/"
    "DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
)

SIGNATURE_FILE = Path(
    "data/signatures/wnt_tcf_signatures.xlsx"
)

OUTPUT_DIR = Path("results/gsea")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------
# Load DESeq2 results
# ---------------------------------------------------------

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

# ---------------------------------------------------------
# Ranking statistic
#
# signed significance score:
# log2FC * -log10(p)
# ---------------------------------------------------------

tiny = np.finfo(float).tiny

df["rank_score"] = (
    df["log2FoldChange"]
    * -np.log10(df["pvalue"].clip(lower=tiny))
)

# If duplicate symbols exist, keep strongest absolute score
ranked = (
    df.sort_values(
        "rank_score",
        key=lambda x: x.abs(),
        ascending=False
    )
    .drop_duplicates("symbol_upper")
    [["symbol_upper", "rank_score"]]
)

# Deterministically order exact ties by gene symbol, then separate them
# at machine precision. This only resolves identical rounded input values.
ranked = ranked.sort_values(["rank_score", "symbol_upper"], ascending=[False, True])
tie_index = ranked.groupby("rank_score", sort=False).cumcount().to_numpy()
ranked["rank_score"] = (
    ranked["rank_score"].to_numpy()
    - tie_index * np.spacing(np.abs(ranked["rank_score"].to_numpy()))
)

ranked = ranked.sort_values(
    "rank_score",
    ascending=False
)

ranked_file = OUTPUT_DIR / "full_transcriptome_ranked.rnk"

ranked.to_csv(
    ranked_file,
    sep="\t",
    index=False,
    header=False
)

print(f"Ranked genes: {len(ranked):,}")
print(f"Saved ranked file: {ranked_file}")

# ---------------------------------------------------------
# Load custom Wnt signatures
# ---------------------------------------------------------

signature_sheets = {
    "BetaCatenin_Independent":
        "bCat indep (65 out of 231)",

    "BetaCatenin_Dependent":
        "bCat dep (166 out of 231)",

    "TCF_Independent":
        "TCF indep (27 out of 166)",

    "TCF_Dependent":
        "TCF dep (139 out of 166)",
}

gene_sets = {}

for name, sheet in signature_sheets.items():

    sig = pd.read_excel(
        SIGNATURE_FILE,
        sheet_name=sheet
    )

    sig.columns = sig.columns.str.strip()

    genes = (
        sig["symbol"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .tolist()
    )

    gene_sets[name] = genes

    print(
        f"{name}: {len(genes)} genes loaded"
    )

# ---------------------------------------------------------
# Add curated PCP set
# ---------------------------------------------------------

gene_sets["WNT_PCP_CORE"] = [
    "WNT5A",
    "ROR2",
    "CELSR1",
    "VANGL2",
    "RAC1",
    "PRICKLE1",
    "FZD7",
    "DVL2",
    "RHOA",
    "ROCK1",
    "MAPK8",
]

# ---------------------------------------------------------
# Add epithelial/junction set
# ---------------------------------------------------------

gene_sets["EPITHELIAL_JUNCTION_CORE"] = [
    "CDH1",
    "EPCAM",
    "OCLN",
    "CLDN3",
    "CLDN4",
    "CLDN6",
    "CLDN7",
    "CLDN8",
    "CLDN9",
    "CRB3",
    "DSP",
    "PKP1",
    "DSG2",
    "DSG3",
    "DSC2",
    "KRT5",
    "KRT14",
    "KRT17",
]

# ---------------------------------------------------------
# Run preranked GSEA
# ---------------------------------------------------------

pre_res = gp.prerank(
    rnk=ranked,
    gene_sets=gene_sets,
    min_size=5,
    max_size=500,
    permutation_num=1000,
        weight=1.0,
        threads=1,
    seed=42,
    outdir=None,
    verbose=True
)

results = pre_res.res2d.copy()

output_file = OUTPUT_DIR / "custom_wnt_preranked_gsea.csv"
results.to_csv(
    output_file,
    index=False
)

print("\n===== PRERANKED GSEA RESULTS =====\n")

show_cols = [
    c for c in [
        "Term",
        "ES",
        "NES",
        "NOM p-val",
        "FDR q-val",
        "Lead_genes"
    ]
    if c in results.columns
]

print(
    results[
        show_cols
    ].sort_values(
        "FDR q-val"
    ).to_string(index=False)
)

print(f"\nSaved: {output_file}")
