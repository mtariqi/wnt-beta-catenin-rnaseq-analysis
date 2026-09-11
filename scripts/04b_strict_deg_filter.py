import pandas as pd
from pathlib import Path

INPUT = Path(
    "data/processed/"
    "DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
)

OUTPUT_DIR = Path("results/deg")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

P_CUTOFF = 0.01
LOG2FC_CUTOFF = 2.0

df = pd.read_csv(INPUT)
df.columns = df.columns.str.strip()

print(f"Raw rows: {len(df):,}")

df = df.dropna(
    subset=["symbol", "log2FoldChange", "pvalue"]
).copy()

strict = df[
    (df["pvalue"] < P_CUTOFF)
    & (df["log2FoldChange"].abs() > LOG2FC_CUTOFF)
].copy()

strict["direction"] = strict["log2FoldChange"].apply(
    lambda x: "Up" if x > 0 else "Down"
)

strict = strict.sort_values(
    ["pvalue", "log2FoldChange"],
    ascending=[True, False]
)

output = OUTPUT_DIR / "strict_DEGs_p0.01_abslog2fc_gt2.csv"
strict.to_csv(output, index=False)

print("\n===== STRICT DEG AUDIT =====")
print(f"Threshold: p < {P_CUTOFF}")
print(f"Threshold: |log2FC| > {LOG2FC_CUTOFF}")
print(f"Strict DEGs: {len(strict):,}")
print(f"Up: {(strict['direction'] == 'Up').sum():,}")
print(f"Down: {(strict['direction'] == 'Down').sum():,}")

print("\nTop genes:")
print(
    strict[
        ["symbol", "log2FoldChange", "pvalue", "direction"]
    ].head(20).to_string(index=False)
)

print(f"\nSaved: {output}")
