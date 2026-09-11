import pandas as pd
from pathlib import Path

INPUT = Path(
    "results/wnt_pcp_mechanistic/"
    "significant_mechanistic_pathways_FDR005.csv"
)

OUTDIR = Path("results/mechanistic_drivers")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# ---------------------------------------------------------
# Clean FDR values
# ---------------------------------------------------------

df["FDR_display"] = df["FDR q-val"].apply(
    lambda x: "<0.001" if x < 0.001 else f"{x:.3f}"
)

# ---------------------------------------------------------
# Direction
# ---------------------------------------------------------

df["Direction"] = df["NES"].apply(
    lambda x: "Positive" if x > 0 else "Negative"
)

# ---------------------------------------------------------
# Select useful columns
# ---------------------------------------------------------

columns = [
    "Mechanistic_Group",
    "Library",
    "Term",
    "NES",
    "FDR q-val",
    "FDR_display",
    "Direction",
]

if "Lead_genes" in df.columns:
    columns.append("Lead_genes")

out = df[columns].copy()

# strongest signals first
out["abs_NES"] = out["NES"].abs()

out = out.sort_values(
    ["Mechanistic_Group", "FDR q-val", "abs_NES"],
    ascending=[True, True, False]
)

out = out.drop(columns="abs_NES")

outfile = OUTDIR / "significant_mechanistic_pathways.csv"

out.to_csv(outfile, index=False)

print("\n===== SIGNIFICANT MECHANISTIC PATHWAYS =====\n")

print(
    out[
        [
            "Mechanistic_Group",
            "Term",
            "NES",
            "FDR_display",
            "Direction",
        ]
    ].to_string(index=False)
)

print(f"\nSaved: {outfile}")

# ---------------------------------------------------------
# Top leading-edge genes
# ---------------------------------------------------------

if "Lead_genes" in df.columns:

    rows = []

    for _, row in df.iterrows():

        genes = str(row["Lead_genes"]).split(";")

        for gene in genes:

            gene = gene.strip()

            if gene and gene.lower() != "nan":

                rows.append(
                    {
                        "Gene": gene,
                        "Mechanistic_Group":
                            row["Mechanistic_Group"],
                        "Term":
                            row["Term"],
                        "NES":
                            row["NES"],
                        "FDR":
                            row["FDR q-val"],
                    }
                )

    genes = pd.DataFrame(
        rows,
        columns=[
            "Gene",
            "Mechanistic_Group",
            "Term",
            "NES",
            "FDR",
        ],
    )

    frequency = (
        genes
        .groupby("Gene")
        .agg(
            Pathway_Count=("Term", "nunique"),
            Mechanistic_Group_Count=(
                "Mechanistic_Group",
                "nunique"
            ),
            Mean_NES=("NES", "mean"),
            Min_FDR=("FDR", "min"),
        )
        .reset_index()
        .sort_values(
            ["Pathway_Count", "Min_FDR"],
            ascending=[False, True]
        )
    )

    gene_file = (
        OUTDIR /
        "recurrent_leading_edge_genes.csv"
    )

    frequency.to_csv(
        gene_file,
        index=False
    )

    print(
        "\n===== RECURRENT LEADING-EDGE GENES =====\n"
    )

    print(
        frequency.head(30).to_string(index=False)
    )

    print(f"\nSaved: {gene_file}")
