#!/usr/bin/env python3
"""Verify human-to-mouse orthology for Doumpas Wnt signatures.

Primary acceptance rule
-----------------------
A mapping is accepted only when:
  * Ensembl Compara reports exactly one mouse target with type
    ``ortholog_one2one``; and
  * the current MGI vertebrate homology report places exactly one human gene
    and exactly one mouse gene in the same homology class; and
  * the Ensembl and MGI mouse symbols agree (case-insensitive).

All other relationships are retained in an audit table with an explicit
exclusion reason. Symbol equality alone is never treated as orthology proof.

The accepted output has ``human_symbol`` and ``mouse_symbol`` columns and can
be supplied directly to ``21_wnt_signature_ranked_validation.py`` using
``--ortholog-map``.

Authoritative sources
---------------------
Ensembl REST/Compara: https://rest.ensembl.org/
MGI report: https://www.informatics.jax.org/downloads/reports/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ENSEMBL_BASE = "https://rest.ensembl.org"
MGI_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
SHEETS = {
    "beta-catenin independent": "bCat indep (65 out of 231)",
    "beta-catenin dependent": "bCat dep (166 out of 231)",
    "TCF/LEF independent": "TCF indep (27 out of 166)",
    "TCF/LEF dependent": "TCF dep (139 out of 166)",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--signatures", required=True, type=Path)
    p.add_argument("--outdir", type=Path,
                   default=Path("results/wnt_orthology_verification"))
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--timeout", type=int, default=60)
    p.add_argument("--retries", type=int, default=4)
    p.add_argument("--refresh", action="store_true",
                   help="Ignore cached Ensembl/MGI responses")
    return p.parse_args()


def http_bytes(url: str, timeout: int, retries: int) -> bytes:
    headers = {"User-Agent": "WntOrthologyAudit/1.0 (academic research)",
               "Accept": "application/json,text/plain,*/*"}
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            # A missing/obsolete symbol is an auditable no-result, not a
            # transient server failure. Returning an empty Ensembl payload
            # prevents repeated requests while keeping the gene excluded.
            if exc.code in {400, 404} and "rest.ensembl.org" in url:
                return b'{"data":[]}'
            last_error = exc
            time.sleep(min(8, 1.5 * (2 ** attempt)))
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
            time.sleep(min(8, 1.5 * (2 ** attempt)))
    raise RuntimeError(f"Request failed after {retries} attempts: {url}: {last_error}")


def load_signatures(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    sheets = set(pd.ExcelFile(path).sheet_names)
    missing = set(SHEETS.values()) - sheets
    if missing:
        raise ValueError(f"Missing signature sheets: {sorted(missing)}")
    rows = []
    for signature, sheet in SHEETS.items():
        tab = pd.read_excel(path, sheet_name=sheet)
        if "symbol" not in tab.columns:
            raise ValueError(f"Sheet {sheet!r} lacks symbol column")
        genes = (tab["symbol"].dropna().astype(str).str.strip().str.upper())
        rows.extend({"signature": signature, "human_symbol": g}
                    for g in genes if g)
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def download_mgi(cache: Path, timeout: int, retries: int, refresh: bool) -> tuple[pd.DataFrame, str]:
    if refresh or not cache.exists():
        cache.write_bytes(http_bytes(MGI_URL, timeout, retries))
    digest = hashlib.sha256(cache.read_bytes()).hexdigest()
    mgi = pd.read_csv(cache, sep="\t", dtype=str, keep_default_na=False)
    expected = {"DB Class Key", "Common Organism Name", "NCBI Taxon ID",
                "Symbol", "Mouse MGI ID", "HGNC ID", "EntrezGene ID"}
    if not expected.issubset(mgi.columns):
        raise ValueError(f"Unexpected MGI columns; missing {sorted(expected-set(mgi.columns))}")
    return mgi, digest


def build_mgi_pairs(mgi: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for class_key, group in mgi.groupby("DB Class Key", sort=False):
        human = group[group["NCBI Taxon ID"].eq("9606")]
        mouse = group[group["NCBI Taxon ID"].eq("10090")]
        n_h, n_m = len(human), len(mouse)
        for _, h in human.iterrows():
            if n_m:
                for _, m in mouse.iterrows():
                    rows.append({
                        "human_symbol": h["Symbol"].upper(),
                        "human_entrez_id": h["EntrezGene ID"],
                        "human_hgnc_id": h["HGNC ID"],
                        "mgi_class_key": class_key,
                        "mgi_human_count": n_h,
                        "mgi_mouse_count": n_m,
                        "mgi_mouse_symbol": m["Symbol"].upper(),
                        "mouse_mgi_id": m["Mouse MGI ID"],
                        "mouse_entrez_id": m["EntrezGene ID"],
                        "mgi_one_to_one": n_h == 1 and n_m == 1,
                    })
            else:
                rows.append({
                    "human_symbol": h["Symbol"].upper(),
                    "human_entrez_id": h["EntrezGene ID"],
                    "human_hgnc_id": h["HGNC ID"],
                    "mgi_class_key": class_key,
                    "mgi_human_count": n_h, "mgi_mouse_count": 0,
                    "mgi_mouse_symbol": "", "mouse_mgi_id": "",
                    "mouse_entrez_id": "", "mgi_one_to_one": False,
                })
    return pd.DataFrame(rows)


def ensembl_homology_url(symbol: str) -> str:
    encoded = urllib.parse.quote(symbol, safe="")
    return (f"{ENSEMBL_BASE}/homology/symbol/homo_sapiens/{encoded}"
            "?target_species=mus_musculus&type=orthologues&aligned=0"
            ";content-type=application/json")


def ensembl_lookup_url(gene_id: str) -> str:
    encoded = urllib.parse.quote(gene_id, safe="")
    return f"{ENSEMBL_BASE}/lookup/id/{encoded}?content-type=application/json"


def cached_json(url: str, cache: Path, timeout: int, retries: int,
                refresh: bool) -> dict:
    if refresh or not cache.exists():
        payload = http_bytes(url, timeout, retries)
        cache.write_bytes(payload)
    return json.loads(cache.read_text(encoding="utf-8"))


def fetch_symbol(symbol: str, cache_dir: Path, timeout: int, retries: int,
                 refresh: bool) -> tuple[str, list[dict], str]:
    cache = cache_dir / f"homology_{urllib.parse.quote(symbol, safe='')}.json"
    try:
        payload = cached_json(ensembl_homology_url(symbol), cache, timeout,
                              retries, refresh)
        data = payload.get("data", [])
        homologies = data[0].get("homologies", []) if data else []
        return symbol, homologies, ""
    except Exception as exc:  # retained in audit rather than silently dropped
        return symbol, [], str(exc)


def lookup_mouse_ids(ids: set[str], cache_dir: Path, workers: int,
                     timeout: int, retries: int, refresh: bool) -> dict[str, str]:
    def one(gene_id: str) -> tuple[str, str]:
        cache = cache_dir / f"lookup_{gene_id}.json"
        try:
            obj = cached_json(ensembl_lookup_url(gene_id), cache, timeout,
                              retries, refresh)
            return gene_id, str(obj.get("display_name", "")).upper()
        except Exception:
            return gene_id, ""
    output = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(one, gene_id) for gene_id in sorted(ids)]
        for future in as_completed(futures):
            key, value = future.result()
            output[key] = value
    return output


def fetch_ensembl(symbols: list[str], cache_dir: Path, workers: int,
                  timeout: int, retries: int, refresh: bool) -> pd.DataFrame:
    results = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fetch_symbol, s, cache_dir, timeout, retries, refresh)
                   for s in symbols]
        for i, future in enumerate(as_completed(futures), 1):
            symbol, homologies, error = future.result()
            if i % 25 == 0 or i == len(futures):
                print(f"Ensembl homology queries: {i}/{len(futures)}")
            if error or not homologies:
                results.append({"human_symbol": symbol, "ensembl_error": error,
                                "ensembl_homology_type": "",
                                "human_ensembl_id": "", "mouse_ensembl_id": "",
                                "human_perc_id": None, "mouse_perc_id": None,
                                "ensembl_taxonomy_level": ""})
                continue
            for h in homologies:
                results.append({
                    "human_symbol": symbol, "ensembl_error": "",
                    "ensembl_homology_type": h.get("type", ""),
                    "human_ensembl_id": h.get("source", {}).get("id", ""),
                    "mouse_ensembl_id": h.get("target", {}).get("id", ""),
                    "human_perc_id": h.get("source", {}).get("perc_id"),
                    "mouse_perc_id": h.get("target", {}).get("perc_id"),
                    "ensembl_taxonomy_level": h.get("taxonomy_level", ""),
                })
    out = pd.DataFrame(results)
    target_ids = set(out["mouse_ensembl_id"].dropna()) - {""}
    names = lookup_mouse_ids(target_ids, cache_dir, workers, timeout, retries, refresh)
    out["ensembl_mouse_symbol"] = out["mouse_ensembl_id"].map(names).fillna("")
    return out


def summarize_ensembl(raw: pd.DataFrame, symbols: list[str]) -> pd.DataFrame:
    rows = []
    for symbol in symbols:
        sub = raw[raw["human_symbol"].eq(symbol)]
        valid = sub[sub["mouse_ensembl_id"].fillna("").ne("")]
        types = sorted(set(valid["ensembl_homology_type"].dropna()) - {""})
        targets = sorted(set(valid["ensembl_mouse_symbol"].dropna()) - {""})
        one = valid[valid["ensembl_homology_type"].eq("ortholog_one2one")]
        rows.append({
            "human_symbol": symbol,
            "ensembl_target_count": len(set(valid["mouse_ensembl_id"])),
            "ensembl_types": ";".join(types),
            "ensembl_mouse_symbols": ";".join(targets),
            "ensembl_one2one": len(valid) == 1 and len(one) == 1,
            "ensembl_mouse_symbol": one.iloc[0]["ensembl_mouse_symbol"] if len(one) == 1 else "",
            "human_ensembl_id": one.iloc[0]["human_ensembl_id"] if len(one) == 1 else "",
            "mouse_ensembl_id": one.iloc[0]["mouse_ensembl_id"] if len(one) == 1 else "",
            "human_perc_id": one.iloc[0]["human_perc_id"] if len(one) == 1 else None,
            "mouse_perc_id": one.iloc[0]["mouse_perc_id"] if len(one) == 1 else None,
            "ensembl_error": ";".join(sorted(set(sub["ensembl_error"].dropna()) - {""})),
        })
    return pd.DataFrame(rows)


def decide(primary_genes: pd.DataFrame, ensembl: pd.DataFrame,
           mgi_pairs: pd.DataFrame) -> pd.DataFrame:
    unique = primary_genes[["human_symbol"]].drop_duplicates()
    audit = unique.merge(ensembl, on="human_symbol", how="left")
    # Collapse MGI to one record only when its class is truly one-to-one.
    mgi_1to1 = (mgi_pairs[mgi_pairs["mgi_one_to_one"]]
                .drop_duplicates("human_symbol", keep=False))
    cols = ["human_symbol", "human_entrez_id", "human_hgnc_id",
            "mgi_class_key", "mgi_human_count", "mgi_mouse_count",
            "mgi_mouse_symbol", "mouse_mgi_id", "mouse_entrez_id",
            "mgi_one_to_one"]
    audit = audit.merge(mgi_1to1[cols], on="human_symbol", how="left")
    audit["mgi_one_to_one"] = audit["mgi_one_to_one"].eq(True)
    audit["database_agreement"] = (
        audit["ensembl_mouse_symbol"].fillna("").str.upper()
        == audit["mgi_mouse_symbol"].fillna("").str.upper()
    ) & audit["ensembl_mouse_symbol"].fillna("").ne("")
    audit["accepted_primary"] = (audit["ensembl_one2one"].fillna(False)
                                 & audit["mgi_one_to_one"]
                                 & audit["database_agreement"])

    def reason(row: pd.Series) -> str:
        if row["accepted_primary"]:
            return "accepted: Ensembl one-to-one and MGI one-to-one agree"
        if row.get("ensembl_error", ""):
            return "excluded: Ensembl query error"
        if not bool(row.get("ensembl_one2one", False)):
            count = int(row.get("ensembl_target_count", 0) or 0)
            return ("excluded: no Ensembl mouse ortholog" if count == 0 else
                    "excluded: Ensembl relationship is not uniquely one-to-one")
        if not bool(row.get("mgi_one_to_one", False)):
            return "excluded: no unique one-human/one-mouse MGI class"
        return "excluded: Ensembl and MGI mouse assignments disagree"

    audit["mapping_decision"] = audit.apply(reason, axis=1)
    audit["mouse_symbol"] = audit["ensembl_mouse_symbol"].where(
        audit["accepted_primary"], "")
    audit = primary_genes.merge(audit, on="human_symbol", how="left")
    return audit.sort_values(["signature", "human_symbol"])


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.workers > 12:
        raise ValueError("--workers must be between 1 and 12")
    args.outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.outdir / "cache"
    cache_dir.mkdir(exist_ok=True)

    signature_rows = load_signatures(args.signatures)
    symbols = sorted(signature_rows["human_symbol"].unique())
    print(f"Unique human signature genes: {len(symbols)}")

    mgi, mgi_sha = download_mgi(cache_dir / "HOM_MouseHumanSequence.rpt",
                                args.timeout, args.retries, args.refresh)
    mgi_pairs = build_mgi_pairs(mgi)
    raw = fetch_ensembl(symbols, cache_dir, args.workers, args.timeout,
                        args.retries, args.refresh)
    ensembl = summarize_ensembl(raw, symbols)
    audit = decide(signature_rows, ensembl, mgi_pairs)

    raw.to_csv(args.outdir / "ensembl_all_candidate_orthologs.csv", index=False)
    audit.to_csv(args.outdir / "wnt_signature_orthology_audit.csv", index=False)
    accepted = (audit[audit["accepted_primary"]]
                [["human_symbol", "mouse_symbol"]].drop_duplicates())
    accepted.to_csv(args.outdir / "human_mouse_consensus_one2one.csv", index=False)
    excluded = audit[~audit["accepted_primary"]]
    excluded.to_csv(args.outdir / "excluded_or_ambiguous_mappings.csv", index=False)

    per_signature = (audit.groupby("signature", as_index=False)
                     .agg(paper_genes=("human_symbol", "nunique"),
                          accepted_consensus_one2one=("accepted_primary", "sum")))
    per_signature["retention_fraction"] = (
        per_signature["accepted_consensus_one2one"] / per_signature["paper_genes"])
    per_signature.to_csv(args.outdir / "orthology_retention_by_signature.csv", index=False)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "ensembl_rest_base": ENSEMBL_BASE,
        "mgi_report_url": MGI_URL,
        "mgi_report_sha256": mgi_sha,
        "acceptance_rule": ("exactly one Ensembl ortholog_one2one target; "
                            "exactly one human and one mouse in MGI class; "
                            "mouse symbols agree"),
        "signature_workbook": str(args.signatures),
        "unique_human_genes": len(symbols),
        "accepted_unique_pairs": len(accepted),
        "excluded_unique_human_genes": len(set(symbols) - set(accepted["human_symbol"])),
    }
    (args.outdir / "orthology_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8")

    print(per_signature.to_string(index=False))
    print(f"Accepted consensus one-to-one pairs: {len(accepted)}/{len(symbols)}")
    print(f"Primary map: {args.outdir/'human_mouse_consensus_one2one.csv'}")
    print(f"Full audit: {args.outdir/'wnt_signature_orthology_audit.csv'}")


if __name__ == "__main__":
    main()
