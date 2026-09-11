import pandas as pd
from pathlib import Path

INPUT = Path(
    "results/pathway_gsea/combined_pathway_preranked_gsea.csv"
)

OUTPUT_DIR = Path(
    "results/wnt_pcp_mechanistic"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# ---------------------------------------------------------
# Normalize text
# ---------------------------------------------------------

df["TERM_UPPER"] = (
    df["Term"]
    .astype(str)
    .str.upper()
)

# ---------------------------------------------------------
# Define mechanistic keyword groups
# ---------------------------------------------------------

keyword_groups = {

    "WNT_BETA_CATENIN": [
        r"\bWNT\b",
        r"WNT SIGNAL",
        r"BETA.?CATENIN",
        r"TCF.?LEF",
        r"CANONICAL WNT",
    ],

    "PCP_POLARITY": [
        r"PLANAR CELL POLAR",
        r"PLANAR POLAR",
        r"PLANAR CELL POLARITY",
        r"NON.?CANONICAL WNT",
        r"WNT.?PCP",
        r"CELL POLARITY",
        r"ESTABLISHMENT OF.*POLARITY",
        r"POLARITY OF",
    ],

    "JUNCTION_ADHESION": [
        r"CELL.?CELL JUNCTION",
        r"TIGHT JUNCTION",
        r"APICAL JUNCTION",
        r"ADHERENS JUNCTION",
        r"DESMOSOM",
        r"CELL.?CELL ADHESION",
        r"CADHERIN",
    ],

    "EPITHELIAL": [
        r"EPITHEL",
        r"KERATIN",
        r"CORNIFIED ENVELOPE",
        r"INTERMEDIATE FILAMENT",
    ],

    "CYTOSKELETON_MIGRATION": [
        r"ACTIN CYTOSKELETON",
        r"CYTOSKELETAL ORGANIZATION",
        r"CELL MIGRATION",
        r"CELL MOTILITY",
        r"RHO GTPASE",
        r"RAC1",
        r"\bPAK",
    ],

    "STRESS_IFN": [
        r"INTERFERON",
        r"CYTOKINE",
        r"INFLAMMAT",
        r"OXIDATIVE STRESS",
        r"CELLULAR STRESS",
        r"\bFOXO",
    ],
}
# ---------------------------------------------------------
# Collect all targeted pathways
# ---------------------------------------------------------

targeted_results = []

for group_name, keywords in keyword_groups.items():

    pattern = "(?:" + "|".join(keywords) + ")"
    subset = df[
        df["TERM_UPPER"]
        .str.contains(
            pattern,
            regex=True,
            na=False
        )
    ].copy()

    subset["Mechanistic_Group"] = group_name

    targeted_results.append(subset)

targeted = pd.concat(
    targeted_results,
    ignore_index=True
)

# remove duplicates that matched multiple keywords
targeted = targeted.drop_duplicates(
    subset=[
        "Library",
        "Term",
        "Mechanistic_Group"
    ]
)

# ---------------------------------------------------------
# Save all targeted results
# ---------------------------------------------------------

targeted = targeted.sort_values(
    [
        "Mechanistic_Group",
        "FDR q-val",
        "NES"
    ]
)

targeted_file = (
    OUTPUT_DIR /
    "all_targeted_mechanistic_pathways.csv"
)

targeted.to_csv(
    targeted_file,
    index=False
)

# ---------------------------------------------------------
# Significant subset
# ---------------------------------------------------------

sig = targeted[
    targeted["FDR q-val"] < 0.05
].copy()

sig_file = (
    OUTPUT_DIR /
    "significant_mechanistic_pathways_FDR005.csv"
)

sig.to_csv(
    sig_file,
    index=False
)

# ---------------------------------------------------------
# Summary by mechanistic group
# ---------------------------------------------------------

summary = (
    targeted
    .groupby("Mechanistic_Group")
    .agg(
        total_pathways=("Term", "count"),
        significant_pathways=(
            "FDR q-val",
            lambda x: (x < 0.05).sum()
        ),
        min_FDR=("FDR q-val", "min"),
        strongest_negative_NES=("NES", "min"),
        strongest_positive_NES=("NES", "max"),
    )
    .reset_index()
)

summary_file = (
    OUTPUT_DIR /
    "mechanistic_group_summary.csv"
)

summary.to_csv(
    summary_file,
    index=False
)

# ---------------------------------------------------------
# Extract leading-edge genes from significant pathways
# ---------------------------------------------------------

leading_edge_rows = []

for _, row in sig.iterrows():

    lead = str(
        row.get("Lead_genes", "")
    )

    if not lead or lead == "nan":
        continue

    genes = [
        g.strip()
        for g in lead.split(";")
        if g.strip()
    ]

    for gene in genes:
        leading_edge_rows.append(
            {
                "Mechanistic_Group":
                    row["Mechanistic_Group"],
                "Library":
                    row["Library"],
                "Term":
                    row["Term"],
                "NES":
                    row["NES"],
                "FDR_q":
                    row["FDR q-val"],
                "Gene":
                    gene,
            }
        )

leading_edges = pd.DataFrame(
    leading_edge_rows,
    columns=[
        "Mechanistic_Group",
        "Library",
        "Term",
        "NES",
        "FDR_q",
        "Gene",
    ],
)

leading_file = (
    OUTPUT_DIR /
    "significant_leading_edge_genes.csv"
)

leading_edges.to_csv(
    leading_file,
    index=False
)

# ---------------------------------------------------------
# Frequency of leading-edge genes
# ---------------------------------------------------------

if not leading_edges.empty:

    gene_frequency = (
        leading_edges
        .groupby(
            [
                "Mechanistic_Group",
                "Gene"
            ]
        )
        .size()
        .reset_index(
            name="Pathway_Count"
        )
        .sort_values(
            [
                "Mechanistic_Group",
                "Pathway_Count"
            ],
            ascending=[
                True,
                False
            ]
        )
    )

else:
    gene_frequency = pd.DataFrame(
        columns=[
            "Mechanistic_Group",
            "Gene",
            "Pathway_Count"
        ]
    )

gene_freq_file = (
    OUTPUT_DIR /
    "leading_edge_gene_frequency.csv"
)

gene_frequency.to_csv(
    gene_freq_file,
    index=False
)

# ---------------------------------------------------------
# Exploratory candidate Wnt and PCP genes
# ---------------------------------------------------------

candidate_genes = {
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
    "CTNNB1",
    "AXIN1",
    "AXIN2",
    "LEF1",
    "TCF7L2",
}

candidate_hits = gene_frequency[
    gene_frequency["Gene"]
    .str.upper()
    .isin(candidate_genes)
].copy()

candidate_hits_file = (
    OUTPUT_DIR /
    "candidate_wnt_pcp_leading_edge_hits.csv"
)

candidate_hits.to_csv(
    candidate_hits_file,
    index=False
)

# ---------------------------------------------------------
# Print concise report
# ---------------------------------------------------------

print("\n===== MECHANISTIC GROUP SUMMARY =====\n")
print(
    summary.to_string(
        index=False
    )
)

print("\n===== SIGNIFICANT WNT / PCP / POLARITY PATHWAYS =====\n")

wnt_pcp_sig = sig[
    sig["Mechanistic_Group"].isin(
        [
            "WNT_BETA_CATENIN",
            "PCP_POLARITY",
        ]
    )
]

show_cols = [
    c for c in [
        "Mechanistic_Group",
        "Library",
        "Term",
        "NES",
        "NOM p-val",
        "FDR q-val",
        "Lead_genes",
    ]
    if c in wnt_pcp_sig.columns
]

if len(wnt_pcp_sig) > 0:

    print(
        wnt_pcp_sig[
            show_cols
        ].sort_values(
            "FDR q-val"
        ).to_string(
            index=False
        )
    )

else:
    print(
        "No Wnt/PCP/polarity pathways reached FDR < 0.05."
    )

print("\n===== EXPLORATORY WNT/PCP CANDIDATE LEADING-EDGE HITS =====\n")

if len(candidate_hits) > 0:
    print(
        candidate_hits.to_string(
            index=False
        )
    )
else:
    print(
        "None of the predefined Wnt/PCP candidate genes "
        "were recurrent leading-edge genes in significant pathways."
    )

print("\nSaved files:")
print(targeted_file)
print(sig_file)
print(summary_file)
print(leading_file)
print(gene_freq_file)
print(candidate_hits_file)
