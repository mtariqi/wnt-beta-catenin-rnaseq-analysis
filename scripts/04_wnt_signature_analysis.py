import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

DESEQ_FILE = Path(
    "data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
)

SIGNATURE_FILE = Path(
    "data/signatures/wnt_tcf_signatures.xlsx"
)

OUTPUT_DIR = Path("results/wnt_signatures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

P_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0

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

df["status"] = "NS"

df.loc[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"] >= LOG2FC_CUTOFF),
    "status"
] = "Up"

df.loc[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"] <= -LOG2FC_CUTOFF),
    "status"
] = "Down"

background_genes = set(df["symbol_upper"])

deg_genes = set(
    df.loc[df["status"] != "NS", "symbol_upper"]
)

# ---------------------------------------------------------
# Signature definitions
# ---------------------------------------------------------

signature_sheets = {
    "beta_catenin_independent": "bCat indep (65 out of 231)",
    "beta_catenin_dependent": "bCat dep (166 out of 231)",
    "tcf_independent": "TCF indep (27 out of 166)",
    "tcf_dependent": "TCF dep (139 out of 166)",
}

summary_rows = []

# ---------------------------------------------------------
# Analyze each signature
# ---------------------------------------------------------

for signature_name, sheet_name in signature_sheets.items():

    sig = pd.read_excel(
        SIGNATURE_FILE,
        sheet_name=sheet_name
    )

    sig.columns = sig.columns.str.strip()

    if "symbol" not in sig.columns:
        raise ValueError(
            f"'symbol' column not found in sheet: {sheet_name}. "
            f"Columns found: {sig.columns.tolist()}"
        )

    signature_genes = (
        sig["symbol"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .drop_duplicates()
    )

    signature_set = set(signature_genes)

    detectable = signature_set & background_genes

    overlap = df[
        df["symbol_upper"].isin(detectable)
    ].copy()

    deg_overlap = overlap[
        overlap["status"] != "NS"
    ].copy()

    up_overlap = overlap[
        overlap["status"] == "Up"
    ].copy()

    down_overlap = overlap[
        overlap["status"] == "Down"
    ].copy()

    # -----------------------------------------------------
    # Fisher exact test
    # -----------------------------------------------------

    a = len(set(deg_overlap["symbol_upper"]))
    b = len(detectable) - a
    c = len(deg_genes) - a
    d = len(background_genes) - a - b - c

    odds_ratio, fisher_p = fisher_exact(
        [[a, b], [c, d]],
        alternative="greater"
    )

    # -----------------------------------------------------
    # Save overlap table
    # -----------------------------------------------------

    overlap_output = (
        OUTPUT_DIR / f"{signature_name}_overlap.csv"
    )

    overlap[
        [
            "gene_id",
            "symbol",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "status"
        ]
    ].sort_values(
        ["status", "pvalue"]
    ).to_csv(
        overlap_output,
        index=False
    )

    summary_rows.append(
        {
            "signature": signature_name,
            "published_genes": len(signature_set),
            "detectable_in_dataset": len(detectable),
            "deg_overlap": len(set(deg_overlap["symbol_upper"])),
            "upregulated_overlap": len(set(up_overlap["symbol_upper"])),
            "downregulated_overlap": len(set(down_overlap["symbol_upper"])),
            "odds_ratio": odds_ratio,
            "fisher_pvalue": fisher_p,
        }
    )

# ---------------------------------------------------------
# Save summary
# ---------------------------------------------------------

summary = pd.DataFrame(summary_rows)
summary["fisher_FDR_BH"] = multipletests(summary["fisher_pvalue"], method="fdr_bh")[1]

summary.to_csv(
    OUTPUT_DIR / "wnt_signature_summary.csv",
    index=False
)

print("\n===== WNT / BETA-CATENIN SIGNATURE SUMMARY =====\n")
print(summary.to_string(index=False))

print(
    "\nSaved results to:",
    OUTPUT_DIR
)
