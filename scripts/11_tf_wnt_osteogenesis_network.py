#!/usr/bin/env python3
"""Local TF-to-Wnt/osteogenesis network analysis for the mouse DEG study.

No study gene symbols are transmitted to external services. Public Enrichr GMT
libraries must already be present in --gene-set-dir.
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


TF_LIBS = ["TRRUST_Transcription_Factors_2019", "ChEA_2022"]
PUBLIC_LIBS = TF_LIBS + ["GO_Biological_Process_2023", "WikiPathways_2024_Mouse"]


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--gene-set-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--p-cutoff", type=float, default=0.05)
    p.add_argument("--top-tfs", type=int, default=8)
    p.add_argument("--genes-per-tf", type=int, default=12)
    p.add_argument("--refresh-public-libraries", action="store_true",
                   help="Redownload public GMT reference libraries before analysis")
    return p.parse_args()


def ensure_public_libraries(directory: Path, refresh: bool = False) -> None:
    """Download public reference gene sets only; never upload study genes."""
    directory.mkdir(parents=True, exist_ok=True)
    missing = [name for name in PUBLIC_LIBS if not (directory / f"{name}.gmt").exists()]
    if not refresh and not missing:
        return
    try:
        import gseapy as gp
    except ImportError as exc:
        raise RuntimeError(
            "gseapy is required to download public gene-set libraries. "
            "Install dependencies with: pip install -r requirements.txt"
        ) from exc
    selected = PUBLIC_LIBS if refresh else missing
    for library in selected:
        gene_sets = gp.get_library(name=library, organism="mouse")
        if not gene_sets:
            raise RuntimeError(f"No records downloaded for public library: {library}")
        with (directory / f"{library}.gmt").open("w", encoding="utf-8") as handle:
            for term, genes in gene_sets.items():
                clean = [str(g).strip().upper() for g in genes if str(g).strip()]
                handle.write("\t".join([term, "public-reference"] + clean) + "\n")


def read_gmt(path: Path) -> dict[str, set[str]]:
    ans = {}
    with path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                ans[parts[0]] = {x.strip().upper() for x in parts[2:] if x.strip()}
    return ans


def tf_name(term: str) -> str:
    return re.split(r"[_ (]", term)[0].upper()


def enrich(query: set[str], universe: set[str], sets: dict[str, set[str]], library: str):
    rows = []
    for term, members0 in sets.items():
        members = members0 & universe
        overlap = query & members
        a = len(overlap)
        if not a:
            continue
        b = len(query - members)
        c = len(members - query)
        d = len(universe - query - members)
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"library": library, "term": term, "tf": tf_name(term),
                     "overlap_count": a, "set_size": len(members),
                     "query_size": len(query), "odds_ratio": odds,
                     "p_value": p, "overlap_genes": ";".join(sorted(overlap))})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = multipletests(out.p_value, method="fdr_bh")[1]
    return out.sort_values(["fdr", "p_value", "tf"])


def pathway_members(go: dict[str, set[str]], wp: dict[str, set[str]]):
    wnt_terms, osteo_terms = {}, {}
    for source, sets in (("GO", go), ("WikiPathways", wp)):
        for term, genes in sets.items():
            low = term.lower()
            if "wnt" in low:
                wnt_terms[f"{source}: {term}"] = genes
            if any(k in low for k in ("ossification", "osteoblast", "osteogenesis",
                                      "bone development", "skeletal system development",
                                      "bone mineralization")):
                osteo_terms[f"{source}: {term}"] = genes
    return wnt_terms, osteo_terms


def union(d):
    return set().union(*d.values()) if d else set()


def collapse_tf(enrichment: pd.DataFrame, n: int):
    return (enrichment.sort_values(["fdr", "p_value"])
            .drop_duplicates("tf").head(n))


def build_edges(top: pd.DataFrame, tf_sets: dict[tuple[str, str], set[str]],
                pathway_genes: set[str], deg: pd.DataFrame, max_genes: int):
    rows = []
    deg_idx = deg.set_index("symbol")
    for _, r in top.iterrows():
        candidates = (tf_sets[(r.library, r.term)] & pathway_genes & set(deg.symbol)) - {r.tf}
        ranked = sorted(candidates, key=lambda g: (deg_idx.loc[g, "pvalue"], -abs(deg_idx.loc[g, "log2FoldChange"])))[:max_genes]
        for gene in ranked:
            x = deg_idx.loc[gene]
            rows.append({"tf": r.tf, "gene": gene, "tf_library": r.library,
                         "tf_enrichment_fdr": r.fdr, "log2FoldChange": x.log2FoldChange,
                         "deg_pvalue": x.pvalue})
    return pd.DataFrame(rows)


def draw_spider(edges: pd.DataFrame, title: str, output: Path):
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    if edges.empty:
        ax.text(.5, .5, "No TF–pathway DEG intersections met the selected criteria",
                ha="center", va="center", fontsize=18)
        ax.axis("off")
        fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig); return
    G = nx.Graph()
    tfs = sorted(edges.tf.unique())
    genes = sorted(edges.gene.unique())
    tf_nodes = [f"TF::{x}" for x in tfs]
    gene_nodes = [f"GENE::{x}" for x in genes]
    G.add_nodes_from(tf_nodes, kind="tf")
    for _, r in edges.iterrows():
        G.add_node(f"GENE::{r.gene}", kind="gene", lfc=r.log2FoldChange)
        G.add_edge(f"TF::{r.tf}", f"GENE::{r.gene}")
    # TFs form an inner ring; target genes form an outer ring.
    pos = {}
    for i, n in enumerate(tf_nodes):
        a = 2*np.pi*i/max(len(tfs),1) + np.pi/2
        pos[n] = (0.42*np.cos(a), 0.42*np.sin(a))
    for i, n in enumerate(gene_nodes):
        a = 2*np.pi*i/max(len(genes),1) + np.pi/2
        pos[n] = (1.0*np.cos(a), 1.0*np.sin(a))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#B7C0CC", width=.9, alpha=.65)
    nx.draw_networkx_nodes(G, pos, nodelist=tf_nodes, ax=ax, node_shape="s",
                           node_color="#7030A0", edgecolors="white", linewidths=1.5, node_size=1350)
    up = [g for g in gene_nodes if G.nodes[g]["lfc"] > 0]
    down = [g for g in gene_nodes if G.nodes[g]["lfc"] < 0]
    nx.draw_networkx_nodes(G, pos, nodelist=up, ax=ax, node_color="#2E8B57",
                           edgecolors="white", linewidths=1, node_size=660)
    nx.draw_networkx_nodes(G, pos, nodelist=down, ax=ax, node_color="#4C78A8",
                           edgecolors="white", linewidths=1, node_size=660)
    nx.draw_networkx_labels(G, pos, labels={n:n.split("::",1)[1] for n in tf_nodes}, ax=ax,
                            font_size=9, font_weight="bold", font_color="white")
    nx.draw_networkx_labels(G, pos, labels={n:n.split("::",1)[1].title() for n in gene_nodes}, ax=ax,
                            font_size=7.5, font_color="#172B4D")
    ax.set_title(title, fontsize=21, weight="bold", color="#15345B", pad=18)
    ax.text(.5, .02, "Squares: enriched TF target sets   Green: upregulated DEG   Blue: downregulated DEG",
            transform=ax.transAxes, ha="center", fontsize=10, color="#475467")
    ax.axis("off"); ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output.with_suffix('.svg'), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    a = args(); a.outdir.mkdir(parents=True, exist_ok=True)
    ensure_public_libraries(a.gene_set_dir, a.refresh_public_libraries)
    df = pd.read_csv(a.input)
    df = df.dropna(subset=["symbol", "pvalue", "log2FoldChange"]).copy()
    df.symbol = df.symbol.astype(str).str.strip().str.upper()
    df = df.sort_values("pvalue").drop_duplicates("symbol")
    universe = set(df.symbol)
    sig = df[df.pvalue < a.p_cutoff]
    up = set(sig.loc[sig.log2FoldChange > 0, "symbol"])
    down = set(sig.loc[sig.log2FoldChange < 0, "symbol"])

    all_enr, tf_sets = [], {}
    for lib in TF_LIBS:
        sets = read_gmt(a.gene_set_dir / f"{lib}.gmt")
        for term, genes in sets.items(): tf_sets[(lib, term)] = genes
        for direction, query in (("downregulated", down), ("upregulated", up)):
            e = enrich(query, universe, sets, lib)
            if not e.empty: e["direction"] = direction; all_enr.append(e)
    enr = pd.concat(all_enr, ignore_index=True).sort_values(["direction", "fdr", "p_value"])
    enr.to_csv(a.outdir / "tf_enrichment_all_local.csv", index=False)

    go = read_gmt(a.gene_set_dir / "GO_Biological_Process_2023.gmt")
    wp = read_gmt(a.gene_set_dir / "WikiPathways_2024_Mouse.gmt")
    wnt_terms, osteo_terms = pathway_members(go, wp)
    wnt = union(wnt_terms) & universe
    osteo = union(osteo_terms) & universe

    down_top = collapse_tf(enr[enr.direction.eq("downregulated")], a.top_tfs)
    up_top = collapse_tf(enr[enr.direction.eq("upregulated")], a.top_tfs)
    # Preserve directionality: the Wnt map uses only downregulated nominal
    # DEGs and the osteogenesis map uses only upregulated nominal DEGs.
    down_df = sig[sig.log2FoldChange < 0].copy()
    up_df = sig[sig.log2FoldChange > 0].copy()
    wnt_edges = build_edges(down_top, tf_sets, wnt, down_df, a.genes_per_tf)
    osteo_edges = build_edges(up_top, tf_sets, osteo, up_df, a.genes_per_tf)
    wnt_edges.to_csv(a.outdir / "wnt_tf_deg_edges.csv", index=False)
    osteo_edges.to_csv(a.outdir / "osteogenesis_tf_deg_edges.csv", index=False)
    draw_spider(wnt_edges, "Downregulated TF programs intersect Wnt signaling genes",
                a.outdir / "wnt_tf_spider_map.png")
    draw_spider(osteo_edges, "Upregulated TF programs intersect osteogenesis genes",
                a.outdir / "osteogenesis_tf_spider_map.png")

    shared = sorted(set(wnt_edges.tf) & set(osteo_edges.tf)) if not wnt_edges.empty and not osteo_edges.empty else []
    pd.DataFrame({"shared_tf": shared}).to_csv(a.outdir / "shared_wnt_osteogenesis_tfs.csv", index=False)
    summary = pd.DataFrame([
        ["All tested unique genes", len(df)], [f"Nominal p < {a.p_cutoff}", len(sig)],
        ["Upregulated nominal DEGs", len(up)], ["Downregulated nominal DEGs", len(down)],
        ["Wnt annotation genes in universe", len(wnt)], ["Osteogenesis annotation genes in universe", len(osteo)],
        ["Wnt TF–DEG edges plotted", len(wnt_edges)], ["Osteogenesis TF–DEG edges plotted", len(osteo_edges)],
        ["TFs shared between maps", len(shared)]], columns=["metric", "value"])
    summary.to_csv(a.outdir / "network_analysis_summary.csv", index=False)
    with (a.outdir / "README.md").open("w") as f:
        f.write("# TF–Wnt–osteogenesis network analysis\n\n")
        f.write("All study-gene calculations were performed locally. TF edges denote membership in public TRRUST/ChEA target sets and overlap with nominal DEGs; they do not establish direct regulation in this tissue. Pathway membership came from GO Biological Process 2023 and WikiPathways 2024 Mouse.\n\n")
        f.write(summary.to_markdown(index=False) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
