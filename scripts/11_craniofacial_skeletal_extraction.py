#!/usr/bin/env python3
"""
Systematic craniofacial/skeletal gene extraction for the Hsu mutant RNA-seq project.

Purpose
-------
1. Load the supplied DESeq2 summary table.
2. Retrieve a current mouse GO Biological Process library through GSEApy/Enrichr.
3. Select craniofacial/skeletal developmental GO terms using an explicit,
   auditable ontology-driven keyword framework.
4. Build the UNION of genes belonging to those terms.
5. Intersect that gene universe with the experimental DESeq2 table.
6. Compute a Benjamini-Hochberg FDR AUDIT because the supplied table lacks
   DESeq2-native padj.
7. Produce tables for Professor Hsu and a prioritized candidate list that can
   be used to build the final Nature-style volcano plot.

IMPORTANT
---------
The BH_FDR_AUDIT column is NOT a substitute for DESeq2-native padj. Replace
it with the original DESeq2 padj when the full DESeq2 results are available.

Run from the repository root:
    python scripts/11_craniofacial_skeletal_extraction.py

Dependencies:
    pandas numpy statsmodels gseapy

Install if needed:
    python -m pip install pandas numpy statsmodels gseapy
"""

from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
import gseapy as gp
from statsmodels.stats.multitest import multipletests


# =====================================================================
# PROJECT PATHS
# =====================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
)

OUTDIR = PROJECT_ROOT / "results" / "craniofacial_skeletal"
OUTDIR.mkdir(parents=True, exist_ok=True)
GO_LIBRARY = "GO_Biological_Process_2023"


# =====================================================================
# ONTOLOGY-DRIVEN TERM FRAMEWORK
# =====================================================================
#
# These categories are intentionally biological-process focused.
# They are NOT hand-selected gene lists.
#
# The exact GO terms available depend on the current Enrichr GO library.
# We therefore retrieve ALL mouse GO BP gene sets and retain terms whose
# names match the process-level patterns below.
#
# Broad enough to capture Professor Hsu's craniofacial/skeletal interests,
# but specific enough to avoid "development" in general.
# =====================================================================

TERM_PATTERNS = {
    "craniofacial_cranial": [
        r"\bcraniofacial\b",
        r"\bcranial skeletal\b",
        r"\bcranial skeleton\b",
        r"\bcranial suture\b",
        r"\bfacial suture\b",
        r"\bfrontonasal\b",
        r"\bface morphogenesis\b",
        r"\bfacial morphogenesis\b",
        r"\bhead morphogenesis\b",
        r"\bcraniofacial morphogenesis\b",
    ],

    "neural_crest": [
        r"\bneural crest\b",
        r"\bcranial neural crest\b",
    ],

    "palate_tooth_jaw": [
        r"\bpalate development\b",
        r"\bpalate morphogenesis\b",
        r"\bpalatal\b",
        r"\btooth development\b",
        r"\btooth morphogenesis\b",
        r"\bodontogenesis\b",
        r"\bmandible\b",
        r"\bmandibular\b",
        r"\bmaxilla\b",
        r"\bmaxillary\b",
        r"\bjaw\b",
    ],

    "skeletal_bone": [
        r"\bskeletal system development\b",
        r"\bskeletal system morphogenesis\b",
        r"\bembryonic skeletal\b",
        r"\bbone development\b",
        r"\bbone morphogenesis\b",
        r"\bossification\b",
        r"\bosteoblast\b",
        r"\bosteogenesis\b",
        r"\bbone mineralization\b",
    ],

    "cartilage_chondrocyte": [
        r"\bcartilage development\b",
        r"\bcartilage morphogenesis\b",
        r"\bchondrocyte\b",
        r"\bchondrogenesis\b",
        r"\bendochondral bone\b",
        r"\bgrowth plate\b",
    ],

    "bmp_developmental": [
        r"\bbmp signaling pathway\b",
        r"\bbone morphogenetic protein\b",
        r"\bresponse to bmp\b",
        r"\bregulation of bmp\b",
    ],
}


# Terms that can match keywords but are usually outside the intended
# craniofacial/skeletal developmental interpretation.
EXCLUDE_PATTERNS = [
    r"cancer",
    r"tumou?r",
    r"immune",
    r"leukocyte",
    r"lymphocyte",
    r"kidney",
    r"renal",
    r"heart",
    r"cardiac",
    r"lung",
    r"liver",
]


# Canonical Wnt genes are NOT used to define the craniofacial/skeletal
# universe. They are retained only as a reference panel for later plots.
CANONICAL_WNT_REFERENCE = [
    "Axin2",
    "Dkk1",
    "Nkd1",
    "Tcf7",
    "Apcdd1",
    "Sp5",
]


# =====================================================================
# HELPERS
# =====================================================================

def choose_mouse_go_bp_library():
    """Choose the newest-looking GO Biological Process library available."""
    libs = gp.get_library_name(organism="mouse")

    candidates = [
        x for x in libs
        if x.startswith("GO_Biological_Process_")
    ]

    if not candidates:
        raise RuntimeError(
            "No GO_Biological_Process library was returned by Enrichr."
        )

    def year_key(name):
        m = re.search(r"(\d{4})$", name)
        return int(m.group(1)) if m else -1

    return sorted(candidates, key=year_key)[-1]


def normalize_term(term):
    return str(term).strip().upper()


def term_is_excluded(term):
    t = str(term).lower()
    return any(re.search(p, t) for p in EXCLUDE_PATTERNS)


def classify_term(term):
    """Return all developmental domains matched by a GO term."""
    t = str(term).lower()

    matches = []

    for domain, patterns in TERM_PATTERNS.items():
        if any(re.search(p, t) for p in patterns):
            matches.append(domain)

    return matches


# =====================================================================
# LOAD DESEQ2 SUMMARY
# =====================================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"\nInput file not found:\n{INPUT_FILE}\n"
        "\nExpected repository layout:\n"
        "data/processed/DESeq2_results_RECALCULATE_4control_4mutant-20251023_TM.csv"
    )

df = pd.read_csv(INPUT_FILE)
df.columns = df.columns.str.strip()

required = {
    "gene_id",
    "symbol",
    "baseMean",
    "log2FoldChange",
    "pvalue",
}

missing = required.difference(df.columns)

if missing:
    raise ValueError(
        f"Missing columns: {sorted(missing)}\n"
        f"Available columns: {df.columns.tolist()}"
    )

print(f"Raw DESeq2 rows: {len(df):,}")

analysis_df = df.dropna(
    subset=["symbol", "log2FoldChange", "pvalue"]
).copy()

analysis_df["symbol"] = analysis_df["symbol"].astype(str).str.strip()
analysis_df["symbol_upper"] = analysis_df["symbol"].str.upper()

# BH audit across all available nominal p-values
analysis_df["BH_FDR_AUDIT"] = multipletests(
    analysis_df["pvalue"].values,
    method="fdr_bh",
)[1]

# Enrichment and candidate tables are symbol-level. Preserve the smallest
# p-value representative where multiple Ensembl IDs map to one symbol.
analysis_df = (
    analysis_df
    .sort_values(["pvalue", "gene_id"], ascending=[True, True])
    .drop_duplicates("symbol_upper", keep="first")
    .copy()
)
print(f"Unique analyzed gene symbols: {len(analysis_df):,}")

analysis_df["Significance"] = "Not significant"

analysis_df.loc[
    (analysis_df["pvalue"] >= 0.05)
    & (analysis_df["pvalue"] < 0.10),
    "Significance"
] = "0.05 ≤ p < 0.10"

analysis_df.loc[
    (analysis_df["pvalue"] < 0.05)
    & (analysis_df["BH_FDR_AUDIT"] >= 0.05),
    "Significance"
] = "Nominal p < 0.05 only"

analysis_df.loc[
    analysis_df["BH_FDR_AUDIT"] < 0.05,
    "Significance"
] = "BH-FDR < 0.05"


# =====================================================================
# RETRIEVE MOUSE GO BIOLOGICAL PROCESS GENE SETS
# =====================================================================

available_libraries = set(gp.get_library_name(organism="mouse"))
if GO_LIBRARY not in available_libraries:
    raise RuntimeError(f"Required Enrichr library not found: {GO_LIBRARY}")

library_name = GO_LIBRARY
print(f"Using pinned mouse library: {library_name}")

go_sets = gp.get_library(
    name=library_name,
    organism="mouse",
)

print(f"GO BP gene sets retrieved: {len(go_sets):,}")

go_snapshot = pd.DataFrame(
    [(term, gene) for term, genes in go_sets.items() for gene in genes],
    columns=["GO_term", "Gene"],
)
go_snapshot.to_csv(
    OUTDIR / "00_GO_Biological_Process_2023_full_membership.csv",
    index=False,
)



# =====================================================================
# SELECT CRANIOFACIAL / SKELETAL TERMS
# =====================================================================

term_records = []
gene_records = []

for term, genes in go_sets.items():

    domains = classify_term(term)

    if not domains:
        continue

    if term_is_excluded(term):
        continue

    clean_genes = sorted(
        {
            str(g).strip()
            for g in genes
            if str(g).strip()
        }
    )

    for domain in domains:
        term_records.append(
            {
                "Domain": domain,
                "GO_term": term,
                "Gene_set_size": len(clean_genes),
                "Library": library_name,
            }
        )

        for gene in clean_genes:
            gene_records.append(
                {
                    "Gene": gene,
                    "Gene_upper": gene.upper(),
                    "Domain": domain,
                    "GO_term": term,
                    "Library": library_name,
                }
            )


terms_df = pd.DataFrame(term_records).drop_duplicates()
genes_long = pd.DataFrame(gene_records).drop_duplicates()

if terms_df.empty:
    raise RuntimeError(
        "No craniofacial/skeletal GO terms matched the configured patterns."
    )

print(f"Selected GO terms: {len(terms_df):,}")
print(
    f"Unique GO-annotated candidate genes: "
    f"{genes_long['Gene_upper'].nunique():,}"
)


# =====================================================================
# BUILD SYSTEMATIC GENE UNIVERSE
# =====================================================================

gene_annotation = (
    genes_long
    .groupby("Gene_upper")
    .agg(
        GO_symbol=("Gene", "first"),
        Domains=(
            "Domain",
            lambda x: ";".join(sorted(set(x)))
        ),
        GO_terms=(
            "GO_term",
            lambda x: " | ".join(sorted(set(x)))
        ),
        Domain_count=(
            "Domain",
            lambda x: len(set(x))
        ),
        GO_term_count=(
            "GO_term",
            lambda x: len(set(x))
        ),
    )
    .reset_index()
)


# =====================================================================
# INTERSECT WITH THE EXPERIMENTAL TRANSCRIPTOME
# =====================================================================

intersection = analysis_df.merge(
    gene_annotation,
    left_on="symbol_upper",
    right_on="Gene_upper",
    how="inner",
)

intersection["abs_log2FC"] = intersection["log2FoldChange"].abs()

# Priority score is for DISPLAY / prioritization only.
# FDR remains the inferential statistic.
intersection["Exploratory_priority_score"] = (
    -np.log10(
        intersection["BH_FDR_AUDIT"].clip(lower=1e-300)
    )
    + 0.25 * intersection["abs_log2FC"]
    + 0.15 * intersection["GO_term_count"]
)

intersection = intersection.sort_values(
    [
        "BH_FDR_AUDIT",
        "pvalue",
        "abs_log2FC",
        "GO_term_count",
    ],
    ascending=[
        True,
        True,
        False,
        False,
    ],
)


# =====================================================================
# DOMAIN SUMMARY
# =====================================================================

domain_rows = []

for domain in TERM_PATTERNS:

    domain_gene_set = set(
        genes_long.loc[
            genes_long["Domain"] == domain,
            "Gene_upper",
        ]
    )

    sub = analysis_df[
        analysis_df["symbol_upper"].isin(domain_gene_set)
    ].copy()

    domain_rows.append(
        {
            "Domain": domain,
            "GO_terms_selected": int(
                terms_df.loc[
                    terms_df["Domain"] == domain,
                    "GO_term",
                ].nunique()
            ),
            "GO_genes_total": len(domain_gene_set),
            "Genes_detected_in_RNAseq": len(sub),
            "BH_FDR_lt_0.05": int(
                (sub["BH_FDR_AUDIT"] < 0.05).sum()
            ),
            "Nominal_p_lt_0.05_only": int(
                (
                    (sub["pvalue"] < 0.05)
                    & (sub["BH_FDR_AUDIT"] >= 0.05)
                ).sum()
            ),
            "p_0.05_to_0.10": int(
                (
                    (sub["pvalue"] >= 0.05)
                    & (sub["pvalue"] < 0.10)
                ).sum()
            ),
            "Median_log2FC": (
                sub["log2FoldChange"].median()
                if len(sub)
                else np.nan
            ),
        }
    )

domain_summary = pd.DataFrame(domain_rows)


# =====================================================================
# PRIORITIZED LIST FOR PROFESSOR / NATURE VOLCANO
# =====================================================================

# Tier 1: all FDR-audit significant genes in the craniofacial/skeletal universe
tier1 = intersection[
    intersection["BH_FDR_AUDIT"] < 0.05
].copy()

# Tier 2: nominal candidates ranked by p, effect size, and annotation recurrence
tier2 = intersection[
    (intersection["pvalue"] < 0.05)
    & (intersection["BH_FDR_AUDIT"] >= 0.05)
].copy()

# Tier 3: borderline 0.05–0.10 candidates
tier3 = intersection[
    (intersection["pvalue"] >= 0.05)
    & (intersection["pvalue"] < 0.10)
].copy()

tier1 = tier1.sort_values(
    ["BH_FDR_AUDIT", "pvalue", "abs_log2FC"],
    ascending=[True, True, False],
)
tier2 = tier2.sort_values(
    ["Exploratory_priority_score", "pvalue"],
    ascending=[False, True],
)
tier3 = tier3.sort_values(
    ["Exploratory_priority_score", "pvalue"],
    ascending=[False, True],
)

tier1["Priority_tier"] = "Tier 1: BH-FDR < 0.05"
tier2["Priority_tier"] = "Tier 2: nominal p < 0.05"
tier3["Priority_tier"] = "Tier 3: 0.05 ≤ p < 0.10"

priority = pd.concat(
    [
        tier1,
        tier2.head(50),
        tier3.head(25),
    ],
    ignore_index=True,
)

# Give preference to genes annotated to multiple relevant GO processes
priority = priority.sort_values(
    [
        "Priority_tier",
        "Exploratory_priority_score",
        "pvalue",
    ],
    ascending=[
        True,
        False,
        True,
    ],
)


# =====================================================================
# CANONICAL WNT REFERENCE TABLE
# =====================================================================

canonical_wnt = analysis_df[
    analysis_df["symbol_upper"].isin(
        {g.upper() for g in CANONICAL_WNT_REFERENCE}
    )
].copy()

canonical_wnt["Panel"] = "Canonical Wnt reference"

canonical_wnt = canonical_wnt[
    [
        "symbol",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "BH_FDR_AUDIT",
        "Significance",
        "Panel",
    ]
].sort_values("pvalue")


# =====================================================================
# SAVE OUTPUTS
# =====================================================================

terms_df.sort_values(
    ["Domain", "GO_term"]
).to_csv(
    OUTDIR / "01_selected_craniofacial_skeletal_GO_terms.csv",
    index=False,
)

genes_long.sort_values(
    ["Domain", "GO_term", "Gene"]
).to_csv(
    OUTDIR / "02_GO_gene_membership_long.csv",
    index=False,
)

gene_annotation.sort_values(
    ["Domain_count", "GO_term_count"],
    ascending=False,
).to_csv(
    OUTDIR / "03_systematic_craniofacial_skeletal_gene_universe.csv",
    index=False,
)

intersection.to_csv(
    OUTDIR / "04_craniofacial_skeletal_DESeq2_intersection.csv",
    index=False,
)

domain_summary.to_csv(
    OUTDIR / "05_domain_summary.csv",
    index=False,
)

priority.to_csv(
    OUTDIR / "06_prioritized_craniofacial_skeletal_candidates.csv",
    index=False,
)

canonical_wnt.to_csv(
    OUTDIR / "07_canonical_Wnt_reference_genes.csv",
    index=False,
)


# =====================================================================
# CONSOLE REPORT
# =====================================================================

print("\n" + "=" * 72)
print("CRANIOFACIAL / SKELETAL EXTRACTION SUMMARY")
print("=" * 72)

print(
    domain_summary.to_string(
        index=False
    )
)

print("\n===== TIER 1: BH-FDR < 0.05 =====\n")

cols = [
    "symbol",
    "log2FoldChange",
    "pvalue",
    "BH_FDR_AUDIT",
    "Domains",
    "GO_term_count",
]

if tier1.empty:
    print("No GO-defined craniofacial/skeletal genes reached BH-FDR < 0.05.")
else:
    print(
        tier1[cols]
        .head(50)
        .to_string(index=False)
    )

print("\n===== TOP TIER 2: NOMINAL p < 0.05 =====\n")

if tier2.empty:
    print("No nominal candidates.")
else:
    print(
        tier2[cols]
        .head(30)
        .to_string(index=False)
    )

print("\n===== CANONICAL WNT REFERENCE =====\n")

print(
    canonical_wnt.to_string(
        index=False
    )
)

print("\nSaved outputs to:")
print(OUTDIR)

print(
    "\nNEXT STEP:\n"
    "Use 06_prioritized_craniofacial_skeletal_candidates.csv "
    "to select labels for the final Nature-style volcano plot."
)
