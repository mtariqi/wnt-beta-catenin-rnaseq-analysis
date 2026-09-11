import pandas as pd
from pathlib import Path

INPUT = Path("data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv")
OUTPUT_DIR = Path("results/deg")

P_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
df.columns = df.columns.str.strip()

print(f"Raw input rows: {len(df):,}")

print("\nMissing values:")
print(
    df[
        ["symbol", "log2FoldChange", "pvalue"]
    ].isna().sum()
)
required = ["gene_id", "symbol", "baseMean", "log2FoldChange", "pvalue"]

missing = [col for col in required if col not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df = df.dropna(
    subset=["symbol", "log2FoldChange", "pvalue"]
).copy()
print(
    f"\nRows available for DEG classification: {len(df):,}"
)
df["regulation"] = "Not significant"

df.loc[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"] >= LOG2FC_CUTOFF),
    "regulation"
] = "Upregulated"

df.loc[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"] <= -LOG2FC_CUTOFF),
    "regulation"
] = "Downregulated"

up = df[df["regulation"] == "Upregulated"].copy()
down = df[df["regulation"] == "Downregulated"].copy()
degs = df[df["regulation"] != "Not significant"].copy()

up = up.sort_values("pvalue")
down = down.sort_values("pvalue")
degs = degs.sort_values("pvalue")

up.to_csv(OUTPUT_DIR / "upregulated_genes.csv", index=False)
down.to_csv(OUTPUT_DIR / "downregulated_genes.csv", index=False)
degs.to_csv(OUTPUT_DIR / "all_DEGs.csv", index=False)

up[["symbol"]].to_csv(
    OUTPUT_DIR / "upregulated_gene_symbols.txt",
    index=False,
    header=False
)

down[["symbol"]].to_csv(
    OUTPUT_DIR / "downregulated_gene_symbols.txt",
    index=False,
    header=False
)

print("===== DEG SUMMARY =====")
print(f"Total genes analysed: {len(df):,}")
print(f"Upregulated: {len(up):,}")
print(f"Downregulated: {len(down):,}")
print(f"Total DEGs: {len(degs):,}")

print("\nTop 10 upregulated genes:")
print(
    up[
        ["symbol", "log2FoldChange", "pvalue"]
    ].head(10).to_string(index=False)
)

print("\nTop 10 downregulated genes:")
print(
    down[
        ["symbol", "log2FoldChange", "pvalue"]
    ].head(10).to_string(index=False)
)
