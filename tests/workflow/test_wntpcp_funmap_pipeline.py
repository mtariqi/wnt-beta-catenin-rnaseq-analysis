import json
import sys
from pathlib import Path

import pandas as pd

PACKAGE = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(PACKAGE))
import wntpcp_funmap_pipeline as pipeline  # noqa: E402


def test_bh_adjust_known_values():
    adjusted = pipeline.bh_adjust([0.01, 0.04, 0.03, 0.002])
    assert list(adjusted.round(3)) == [0.02, 0.04, 0.04, 0.008]


def test_end_to_end_outputs_and_layer_separation(tmp_path, monkeypatch):
    de = tmp_path / "deseq.csv"
    orth = tmp_path / "orthology.csv"
    funmap = tmp_path / "funmap.tsv"
    string = tmp_path / "string.tsv"
    literature = tmp_path / "literature.csv"
    out = tmp_path / "out"

    pd.DataFrame({
        "gene_id": ["m1", "m2", "m3"],
        "symbol": ["Wnt5a", "Ror2", "Dsp"],
        "baseMean": [100, 80, 120],
        "log2FoldChange": [1.2, 0.8, -2.1],
        "pvalue": [0.01, 0.04, 0.001],
    }).to_csv(de, index=False)
    pd.DataFrame({
        "human_symbol": ["WNT5A", "ROR2", "DSP"],
        "mouse_symbol": ["Wnt5a", "Ror2", "Dsp"],
    }).to_csv(orth, index=False)
    pd.DataFrame({
        "gene1": ["WNT5A", "WNT5A"],
        "gene2": ["ROR2", "DSP"],
        "score": [0.93, 0.61],
    }).to_csv(funmap, sep="\t", index=False)
    pd.DataFrame({
        "preferredName_A": ["Wnt5a"], "preferredName_B": ["Dsp"],
        "combined_score": [850],
    }).to_csv(string, sep="\t", index=False)
    pd.DataFrame({
        "source": ["Ror2"], "target": ["Dsp"], "species": ["mouse"],
        "relationship": ["curated association"], "citation": ["PMID:123"],
    }).to_csv(literature, index=False)

    monkeypatch.setattr(sys, "argv", [
        "wntpcp_funmap_pipeline.py", "--deseq", str(de), "--orthology", str(orth),
        "--funmap", str(funmap), "--string", str(string), "--literature", str(literature),
        "--genes", "Wnt5a,Ror2,Dsp", "--strict-orthology", "--outdir", str(out),
    ])
    pipeline.main()

    for suffix in ["png", "pdf", "svg"]:
        assert (out / f"wntpcp_funmap_network.{suffix}").is_file()
    edges = pd.read_csv(out / "wntpcp_funmap_network_edge_audit.csv")
    assert set(edges["evidence_layers"]) == {"FunMap", "FunMap+STRING", "literature"}
    prov = json.loads((out / "wntpcp_funmap_network_run_provenance.json").read_text())
    assert prov["funmap"]["source_format"] == "official_three_column_funmap"
    assert prov["orthology"]["unmapped_genes"] == []
