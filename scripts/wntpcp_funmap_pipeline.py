#!/usr/bin/env python3
"""Auditable overlay of mouse DESeq2 results on human FunMap associations.

This program deliberately keeps four evidence layers separate:
1. current-study mouse differential expression;
2. consensus one-to-one mouse-human orthology;
3. imported FunMap predicted functional-association scores; and
4. optional STRING and literature evidence.

FunMap scores are imported, never reconstructed from DESeq2 summary statistics.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd


DEFAULT_MOUSE_GENES = [
    "Wnt5a", "Vangl2", "Celsr1", "Ror2", "Rac1", "Prickle1", "Egr3",
    "Dsp", "Cdh1", "Krt14", "Krt5", "Bmp6", "Lama5", "Gm6565", "Il31ra",
]

REQUIRED_DESEQ = {"symbol", "log2FoldChange", "pvalue"}
REQUIRED_ORTHOLOGY = {"human_symbol", "mouse_symbol"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Overlay mouse DESeq2 results on orthology-mapped human FunMap edges."
    )
    p.add_argument("--deseq", required=True, type=Path,
                   help="Authoritative mouse DESeq2 CSV/TSV with symbol, log2FoldChange, pvalue.")
    p.add_argument("--orthology", required=True, type=Path,
                   help="Consensus one-to-one CSV/TSV with human_symbol and mouse_symbol.")
    p.add_argument("--funmap", required=True, type=Path,
                   help="Official FunMap gene1/gene2/score file, or a documented derived table.")
    p.add_argument("--string", dest="string_file", type=Path,
                   help="Optional STRING edge file. See README for accepted columns.")
    p.add_argument("--string-species", choices=["mouse", "human"], default="mouse")
    p.add_argument("--literature", type=Path,
                   help="Optional curated edge table with source,target,citation.")
    p.add_argument("--genes", type=str,
                   help="Comma-separated mouse symbols. Defaults to the Wnt/PCP panel.")
    p.add_argument("--gene-file", type=Path,
                   help="Optional one-column mouse-symbol file; overrides --genes/default.")
    p.add_argument("--funmap-score-min", type=float, default=0.0,
                   help="Minimum imported FunMap probability/score to display (default: 0).")
    p.add_argument("--string-score-min", type=float, default=0.0,
                   help="Minimum normalized STRING confidence in [0,1] (default: 0).")
    p.add_argument("--max-edges", type=int, default=80,
                   help="Maximum displayed edges after deterministic evidence-aware ranking.")
    p.add_argument("--strict-orthology", action="store_true",
                   help="Stop if any requested mouse gene lacks a one-to-one human mapping.")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--prefix", default="wntpcp_funmap_network")
    p.add_argument("--title", default="Wnt/PCP functional-association network")
    p.add_argument("--seed", type=int, default=101)
    return p.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    """Read comma/tab-delimited text (optionally gzipped) with delimiter detection."""
    if not path.is_file():
        raise FileNotFoundError(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
        sample = handle.read(8192)
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in sample else ","
    return pd.read_csv(path, sep=delimiter, compression="infer", low_memory=False)


def require_columns(df: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def upper(value: object) -> str:
    return str(value).strip().upper() if pd.notna(value) else ""


def bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def load_deseq(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = read_table(path)
    require_columns(raw, REQUIRED_DESEQ, "DESeq2 table")
    df = raw.copy()
    df["symbol_upper"] = df["symbol"].map(upper)
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    valid = df["pvalue"].between(0, 1, inclusive="both") & df["symbol_upper"].ne("")
    if not valid.any():
        raise ValueError("DESeq2 table contains no valid symbol/pvalue rows")

    # BH is calculated over the complete valid gene-level testing universe,
    # before any display-oriented symbol deduplication.
    recalculated = pd.Series(np.nan, index=df.index, dtype=float)
    recalculated.loc[valid] = bh_adjust(df.loc[valid, "pvalue"].to_numpy())
    if "padj" in df.columns:
        source_padj = pd.to_numeric(df["padj"], errors="coerce")
        df["padj_plot"] = source_padj.where(source_padj.notna(), recalculated)
        padj_source = "DESeq2 padj where present; otherwise BH recomputed over valid p-values"
    else:
        df["padj_plot"] = recalculated
        padj_source = "BH recomputed over all valid p-values before symbol deduplication"

    df["_abs_lfc"] = df["log2FoldChange"].abs().fillna(-1)
    sort_cols = ["pvalue", "_abs_lfc"]
    ascending = [True, False]
    if "baseMean" in df.columns:
        df["baseMean"] = pd.to_numeric(df["baseMean"], errors="coerce")
        sort_cols.append("baseMean")
        ascending.append(False)
    dedup = (df.loc[valid]
               .sort_values(sort_cols, ascending=ascending, kind="mergesort")
               .drop_duplicates("symbol_upper", keep="first")
               .drop(columns=["_abs_lfc"]))
    meta = {
        "input_rows": int(len(raw)),
        "valid_testing_rows": int(valid.sum()),
        "unique_display_symbols": int(dedup["symbol_upper"].nunique()),
        "padj_plot_definition": padj_source,
    }
    return dedup, meta


def load_orthology(path: Path) -> pd.DataFrame:
    df = read_table(path)
    require_columns(df, REQUIRED_ORTHOLOGY, "Orthology table")
    out = df[["human_symbol", "mouse_symbol"]].dropna().copy()
    out["human_upper"] = out["human_symbol"].map(upper)
    out["mouse_upper"] = out["mouse_symbol"].map(upper)
    out = out[(out.human_upper != "") & (out.mouse_upper != "")].drop_duplicates()
    hdup = out.loc[out.duplicated("human_upper", keep=False), "human_symbol"].unique()
    mdup = out.loc[out.duplicated("mouse_upper", keep=False), "mouse_symbol"].unique()
    if len(hdup) or len(mdup):
        raise ValueError(
            "Orthology input is not one-to-one after normalization; "
            f"duplicated human={list(hdup[:10])}, mouse={list(mdup[:10])}"
        )
    return out


def pick_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    by_lower = {c.lower(): c for c in df.columns}
    return next((by_lower[x.lower()] for x in candidates if x.lower() in by_lower), None)


def load_funmap(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    df = read_table(path)
    a = pick_column(df, ["gene1", "gene_a", "source"])
    b = pick_column(df, ["gene2", "gene_b", "target"])
    score = pick_column(df, ["score", "probability", "funmap_score"])
    if not (a and b and score):
        raise ValueError(
            "FunMap table must contain gene1/gene2/score (official output) or "
            "Gene1/Gene2/Probability (documented derived table)."
        )
    official = a.lower() == "gene1" and b.lower() == "gene2" and score.lower() == "score"
    out = pd.DataFrame({
        "human_a": df[a].map(upper),
        "human_b": df[b].map(upper),
        "funmap_score": pd.to_numeric(df[score], errors="coerce"),
    })
    llr = pick_column(df, ["LLR_Score", "llr", "log_likelihood_ratio"])
    out["funmap_llr"] = pd.to_numeric(df[llr], errors="coerce") if llr else np.nan
    out["funmap_source_format"] = (
        "official_three_column_funmap" if official else "derived_funmap_probability_table"
    )
    out = out.dropna(subset=["funmap_score"])
    out = out[(out.human_a != "") & (out.human_b != "") & (out.human_a != out.human_b)]
    out[["human_a", "human_b"]] = np.sort(out[["human_a", "human_b"]].to_numpy(), axis=1)
    out = (out.sort_values(["funmap_score", "funmap_llr"], ascending=False, kind="mergesort")
              .drop_duplicates(["human_a", "human_b"], keep="first"))
    meta = {
        "source_format": out["funmap_source_format"].iloc[0] if len(out) else "empty",
        "score_column": score,
        "score_interpretation": "imported predicted functional-association score/probability; not a p-value",
        "input_edges": int(len(df)),
        "valid_unique_edges": int(len(out)),
    }
    return out, meta


def canonical_mouse_edges(df: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    out = df.copy()
    out["mouse_a"] = out[a].map(upper)
    out["mouse_b"] = out[b].map(upper)
    out = out[(out.mouse_a != "") & (out.mouse_b != "") & (out.mouse_a != out.mouse_b)]
    out[["mouse_a", "mouse_b"]] = np.sort(out[["mouse_a", "mouse_b"]].to_numpy(), axis=1)
    return out


def load_string(path: Path, species: str, orth: pd.DataFrame) -> pd.DataFrame:
    df = read_table(path)
    a = pick_column(df, ["preferredName_A", "source", "gene1", "protein1"])
    b = pick_column(df, ["preferredName_B", "target", "gene2", "protein2"])
    score = pick_column(df, ["score", "combined_score", "string_score"])
    if not (a and b and score):
        raise ValueError("STRING table needs endpoint columns and score/combined_score")
    tmp = pd.DataFrame({a: df[a], b: df[b], "string_score": pd.to_numeric(df[score], errors="coerce")})
    if tmp["string_score"].max(skipna=True) > 1:
        tmp["string_score"] = tmp["string_score"] / 1000.0
    if species == "human":
        h2m = orth.set_index("human_upper")["mouse_upper"]
        tmp[a] = tmp[a].map(upper).map(h2m)
        tmp[b] = tmp[b].map(upper).map(h2m)
    out = canonical_mouse_edges(tmp, a, b)
    out = out.dropna(subset=["string_score"])
    return (out.sort_values("string_score", ascending=False, kind="mergesort")
               .drop_duplicates(["mouse_a", "mouse_b"], keep="first")
               [["mouse_a", "mouse_b", "string_score"]])


def load_literature(path: Path, orth: pd.DataFrame) -> pd.DataFrame:
    df = read_table(path)
    require_columns(df, {"source", "target", "citation"}, "Literature edge table")
    out = df.copy()
    if "species" not in out.columns:
        out["species"] = "mouse"
    h2m = orth.set_index("human_upper")["mouse_upper"]
    human = out["species"].astype(str).str.lower().eq("human")
    out.loc[human, "source"] = out.loc[human, "source"].map(upper).map(h2m)
    out.loc[human, "target"] = out.loc[human, "target"].map(upper).map(h2m)
    out = canonical_mouse_edges(out, "source", "target")
    for col in ["relationship", "evidence_level", "notes"]:
        if col not in out.columns:
            out[col] = ""
    agg = (out.groupby(["mouse_a", "mouse_b"], as_index=False)
              .agg(literature_citation=("citation", lambda s: " | ".join(sorted(set(map(str, s))))),
                   literature_relationship=("relationship", lambda s: " | ".join(sorted(set(filter(None, map(str, s)))))),
                   literature_evidence_level=("evidence_level", lambda s: " | ".join(sorted(set(filter(None, map(str, s)))))),
                   literature_notes=("notes", lambda s: " | ".join(sorted(set(filter(None, map(str, s))))))))
    return agg


def requested_genes(args: argparse.Namespace) -> list[str]:
    if args.gene_file:
        lines = [x.strip() for x in args.gene_file.read_text(encoding="utf-8").splitlines()]
        genes = [x for x in lines if x and not x.startswith("#")]
    elif args.genes:
        genes = [x.strip() for x in args.genes.split(",") if x.strip()]
    else:
        genes = DEFAULT_MOUSE_GENES
    seen: set[str] = set()
    return [g for g in genes if not (upper(g) in seen or seen.add(upper(g)))]


def build_node_audit(genes: list[str], de: pd.DataFrame, orth: pd.DataFrame) -> pd.DataFrame:
    requested = pd.DataFrame({"requested_mouse_symbol": genes})
    requested["mouse_upper"] = requested.requested_mouse_symbol.map(upper)
    o = orth[["mouse_upper", "mouse_symbol", "human_symbol", "human_upper"]]
    keep = [c for c in ["symbol", "symbol_upper", "gene_id", "baseMean", "log2FoldChange", "pvalue", "padj_plot"] if c in de.columns]
    d = de[keep].rename(columns={"symbol_upper": "mouse_upper", "symbol": "deseq_mouse_symbol"})
    out = requested.merge(o, on="mouse_upper", how="left").merge(d, on="mouse_upper", how="left")
    out["orthology_status"] = np.where(out.human_upper.notna(), "mapped_one_to_one", "unmapped")
    out["deseq_status"] = np.where(out.log2FoldChange.notna(), "represented", "not_represented")
    return out


def build_edge_audit(
    nodes: pd.DataFrame,
    funmap: pd.DataFrame,
    string_edges: pd.DataFrame | None,
    literature: pd.DataFrame | None,
    funmap_min: float,
    string_min: float,
) -> pd.DataFrame:
    h2m = nodes.dropna(subset=["human_upper"]).set_index("human_upper")["mouse_upper"]
    fm = funmap.copy()
    fm["mouse_a"] = fm.human_a.map(h2m)
    fm["mouse_b"] = fm.human_b.map(h2m)
    fm = fm.dropna(subset=["mouse_a", "mouse_b"])
    fm[["mouse_a", "mouse_b"]] = np.sort(fm[["mouse_a", "mouse_b"]].to_numpy(), axis=1)
    fm = fm[fm.funmap_score >= funmap_min]
    fm = fm[["mouse_a", "mouse_b", "human_a", "human_b", "funmap_score", "funmap_llr", "funmap_source_format"]]

    frames = [fm[["mouse_a", "mouse_b"]]]
    if string_edges is not None:
        string_edges = string_edges[string_edges.string_score >= string_min].copy()
        frames.append(string_edges[["mouse_a", "mouse_b"]])
    if literature is not None:
        frames.append(literature[["mouse_a", "mouse_b"]])
    keys = pd.concat(frames, ignore_index=True).drop_duplicates()
    out = keys.merge(fm, on=["mouse_a", "mouse_b"], how="left")
    if string_edges is not None:
        out = out.merge(string_edges, on=["mouse_a", "mouse_b"], how="left")
    else:
        out["string_score"] = np.nan
    if literature is not None:
        out = out.merge(literature, on=["mouse_a", "mouse_b"], how="left")
    else:
        for c in ["literature_citation", "literature_relationship", "literature_evidence_level", "literature_notes"]:
            out[c] = ""

    out["funmap_present"] = out.funmap_score.notna()
    out["string_present"] = out.string_score.notna()
    out["literature_present"] = out.literature_citation.fillna("").ne("")
    codes = []
    for row in out.itertuples():
        layer = []
        if row.funmap_present: layer.append("FunMap")
        if row.string_present: layer.append("STRING")
        if row.literature_present: layer.append("literature")
        codes.append("+".join(layer))
    out["evidence_layers"] = codes
    out["evidence_layer_count"] = out[["funmap_present", "string_present", "literature_present"]].sum(axis=1)

    stats = nodes.set_index("mouse_upper")
    m2h = nodes.dropna(subset=["human_upper"]).set_index("mouse_upper")["human_symbol"]
    for side in ["a", "b"]:
        key = out[f"mouse_{side}"]
        out[f"mouse_symbol_{side}"] = key.map(nodes.set_index("mouse_upper")["requested_mouse_symbol"])
        out[f"human_symbol_{side}"] = key.map(m2h)
        for stat in ["log2FoldChange", "pvalue", "padj_plot"]:
            if stat in stats.columns:
                out[f"gene_{side}_{stat}"] = key.map(stats[stat])
    out["display_edge"] = True
    out["interpretation"] = (
        "Association evidence only; no edge establishes direct binding, directionality, "
        "regulator activity, or causality in the current experiment."
    )
    return out.sort_values(
        ["evidence_layer_count", "funmap_score", "string_score"],
        ascending=[False, False, False], kind="mergesort"
    ).reset_index(drop=True)


def select_edges(edge_audit: pd.DataFrame, max_edges: int) -> pd.DataFrame:
    ranked = edge_audit.copy()
    ranked["_fun"] = ranked.funmap_score.fillna(-1)
    ranked["_string"] = ranked.string_score.fillna(-1)
    ranked = ranked.sort_values(
        ["evidence_layer_count", "_fun", "_string", "mouse_a", "mouse_b"],
        ascending=[False, False, False, True, True], kind="mergesort"
    )
    return ranked.head(max_edges).drop(columns=["_fun", "_string"])


def draw_network(nodes: pd.DataFrame, edges: pd.DataFrame, outdir: Path, prefix: str,
                 title: str, seed: int) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 11, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
    })
    G = nx.Graph()
    node_lookup = nodes.set_index("mouse_upper")
    for row in nodes.itertuples():
        if row.deseq_status == "represented":
            G.add_node(row.mouse_upper)
    for row in edges.itertuples():
        if row.mouse_a in G and row.mouse_b in G:
            G.add_edge(row.mouse_a, row.mouse_b)
    if G.number_of_nodes() == 0:
        raise ValueError("No requested genes are represented in the DESeq2 table")

    pos = nx.spring_layout(G, seed=seed, k=max(0.55, 2.0 / math.sqrt(max(2, G.number_of_nodes()))))
    fig, ax = plt.subplots(figsize=(11.2, 8.2), constrained_layout=True)

    # Evidence-specific edge encoding. The audit CSV remains the authoritative record.
    styles = {
        "FunMap": ("#2563A6", "solid", 1.8),
        "STRING": ("#6B7280", (0, (5, 3)), 1.6),
        "literature": ("#B7791F", (0, (1, 2)), 1.8),
        "multi": ("#6B4C9A", "solid", 2.4),
    }
    for row in edges.itertuples():
        if not G.has_edge(row.mouse_a, row.mouse_b):
            continue
        key = "multi" if row.evidence_layer_count > 1 else row.evidence_layers
        color, style, base_width = styles[key]
        strength = row.funmap_score if pd.notna(row.funmap_score) else row.string_score
        width = base_width + (float(strength) if pd.notna(strength) else 0.0)
        nx.draw_networkx_edges(G, pos, edgelist=[(row.mouse_a, row.mouse_b)],
                               edge_color=color, style=style, width=width, alpha=0.9, ax=ax)

    lfcs = np.array([float(node_lookup.loc[n, "log2FoldChange"]) for n in G.nodes()])
    finite = np.abs(lfcs[np.isfinite(lfcs)])
    vmax = max(1.0, float(np.quantile(finite, 0.95)) if finite.size else 1.0)
    node_collection = nx.draw_networkx_nodes(
        G, pos, node_color=lfcs, cmap="coolwarm", vmin=-vmax, vmax=vmax,
        node_size=1650, edgecolors="#1F2937", linewidths=1.1, ax=ax,
    )
    labels = {n: str(node_lookup.loc[n, "requested_mouse_symbol"]) for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8.5, font_weight="bold", ax=ax)

    # Compact evidence-code labels; full scores and citations are in the edge audit.
    edge_labels = {}
    for row in edges.itertuples():
        if G.has_edge(row.mouse_a, row.mouse_b):
            codes = []
            if row.funmap_present: codes.append("F")
            if row.string_present: codes.append("S")
            if row.literature_present: codes.append("L")
            edge_labels[(row.mouse_a, row.mouse_b)] = "+".join(codes)
    nx.draw_networkx_edge_labels(
        G, pos, edge_labels=edge_labels, font_size=6.2, rotate=False,
        bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.82}, ax=ax,
    )

    cbar = fig.colorbar(node_collection, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(r"Current-study $log_2$ fold change (mutant vs control)")
    legend = [
        Line2D([0], [0], color=styles["FunMap"][0], lw=2, label="F: FunMap association"),
        Line2D([0], [0], color=styles["STRING"][0], lw=2, ls="--", label="S: STRING association"),
        Line2D([0], [0], color=styles["literature"][0], lw=2, ls=":", label="L: curated literature"),
        Line2D([0], [0], color=styles["multi"][0], lw=2.5, label="Multiple evidence layers"),
    ]
    ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=4, frameon=False)
    ax.set_title(title + "\nMouse DESeq2 nodes mapped to human one-to-one orthologs for FunMap lookup")
    ax.text(
        0.5, -0.13,
        "Edges are undirected association evidence. FunMap/STRING scores are not p-values; literature edges are curated. "
        "No edge establishes direct binding, activity, directionality, or causality.",
        transform=ax.transAxes, ha="center", va="top", fontsize=7, color="#4B5563", wrap=True,
    )
    ax.axis("off")
    for ext, kwargs in [("png", {"dpi": 600}), ("pdf", {}), ("svg", {})]:
        fig.savefig(outdir / f"{prefix}.{ext}", bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    genes = requested_genes(args)
    de, de_meta = load_deseq(args.deseq)
    orth = load_orthology(args.orthology)
    funmap, funmap_meta = load_funmap(args.funmap)
    nodes = build_node_audit(genes, de, orth)
    unmapped = nodes.loc[nodes.orthology_status == "unmapped", "requested_mouse_symbol"].tolist()
    if args.strict_orthology and unmapped:
        raise ValueError(f"Requested genes lacking one-to-one human orthologs: {unmapped}")

    string_edges = load_string(args.string_file, args.string_species, orth) if args.string_file else None
    literature = load_literature(args.literature, orth) if args.literature else None
    selected_mouse = set(nodes.mouse_upper)
    if string_edges is not None:
        string_edges = string_edges[string_edges.mouse_a.isin(selected_mouse) & string_edges.mouse_b.isin(selected_mouse)]
    if literature is not None:
        literature = literature[literature.mouse_a.isin(selected_mouse) & literature.mouse_b.isin(selected_mouse)]
    edge_audit = build_edge_audit(
        nodes, funmap, string_edges, literature,
        args.funmap_score_min, args.string_score_min,
    )
    displayed = select_edges(edge_audit, args.max_edges)
    edge_audit["display_edge"] = edge_audit.index.isin(displayed.index)

    nodes.to_csv(args.outdir / f"{args.prefix}_node_audit.csv", index=False)
    edge_audit.to_csv(args.outdir / f"{args.prefix}_edge_audit.csv", index=False)
    draw_network(nodes, displayed, args.outdir, args.prefix, args.title, args.seed)

    inputs = [args.deseq, args.orthology, args.funmap]
    inputs += [p for p in [args.string_file, args.literature] if p]
    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(map(str, sys.argv)),
        "python": platform.python_version(),
        "packages": {"pandas": pd.__version__, "numpy": np.__version__,
                     "matplotlib": mpl.__version__, "networkx": nx.__version__},
        "inputs": [{"path": str(p.resolve()), "sha256": sha256(p)} for p in inputs],
        "parameters": {
            "genes": genes, "funmap_score_min": args.funmap_score_min,
            "string_score_min": args.string_score_min, "max_edges": args.max_edges,
            "strict_orthology": args.strict_orthology, "seed": args.seed,
        },
        "deseq": de_meta,
        "orthology": {
            "mapping_rows": int(len(orth)), "requested": len(genes),
            "mapped": int((nodes.orthology_status == "mapped_one_to_one").sum()),
            "unmapped_genes": unmapped,
        },
        "funmap": funmap_meta,
        "outputs": {
            "candidate_edges": int(len(edge_audit)), "displayed_edges": int(len(displayed)),
            "represented_nodes": int((nodes.deseq_status == "represented").sum()),
        },
        "interpretation_constraints": [
            "FunMap scores are imported functional-association predictions, not gene-level p-values.",
            "STRING scores are association confidence scores, not p-values.",
            "Literature edges are curated claims and remain distinct from FunMap and STRING.",
            "No edge demonstrates direct binding, regulator activity, directionality, or causality.",
            "Current-study tissue specificity applies to node differential expression, not imported FunMap edges.",
        ],
    }
    (args.outdir / f"{args.prefix}_run_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote auditable outputs to {args.outdir}")
    if unmapped:
        print(f"WARNING: {len(unmapped)} requested genes lacked one-to-one mappings; see node audit.")


if __name__ == "__main__":
    main()
