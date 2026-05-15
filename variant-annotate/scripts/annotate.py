#!/usr/bin/env python3
"""variant-annotate: per-rsID lookup against gnomAD v4 + GTEx v8.

Given a list of rsIDs, fetches:
  - gnomAD: most-severe coding consequence, gene, per-pop allele freq
  - GTEx:   top single-tissue eQTL and sQTL (smallest p-value)

Output is a TSV (default), markdown table, or JSON. Per-rsID JSON is
cached under ~/.cache/variant-annotate/ so repeated calls are free.

See SKILL.md for full usage. CLI in __main__ below."""

from __future__ import annotations
import argparse
import json
import sys
import urllib.request
from pathlib import Path
from typing import Optional

DEFAULT_CACHE = Path.home() / ".cache" / "variant-annotate"
GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
GTEX_BASE = "https://gtexportal.org/api/v2"
UA = "statgen-skills/variant-annotate (annotate.py)"

# Population columns we keep from gnomAD (top-level, no sex subpops).
GNOMAD_POPS = ["nfe", "afr", "eas", "sas", "amr", "fin", "asj", "mid",
               "remaining"]

# Most-severe consequence ordering (subset; covers cases that matter
# for fine-mapping interpretation — the lower the rank, the more
# severe). Anything not listed gets rank 99.
CSQ_SEVERITY = {
    "transcript_ablation": 0, "splice_acceptor_variant": 1,
    "splice_donor_variant": 2, "stop_gained": 3,
    "frameshift_variant": 4, "stop_lost": 5, "start_lost": 6,
    "missense_variant": 7, "splice_region_variant": 8,
    "synonymous_variant": 9, "5_prime_UTR_variant": 10,
    "3_prime_UTR_variant": 11, "intron_variant": 12,
    "upstream_gene_variant": 13, "downstream_gene_variant": 14,
    "non_coding_transcript_exon_variant": 15,
    "regulatory_region_variant": 16,
    "intergenic_variant": 17,
}

_GNOMAD_QUERY = """
query ($rsid: String!) {
  variant(rsid: $rsid, dataset: gnomad_r4) {
    variant_id
    rsids
    transcript_consequences {
      gene_symbol
      consequence_terms
      hgvsc
      hgvsp
    }
    exome { populations { id ac an } }
    genome { populations { id ac an } }
  }
}
"""


def gnomad_query(rsid: str, cache_dir: Path,
                  use_cache: bool = True) -> Optional[dict]:
    """Fetch gnomAD v4 annotation for an rsID. Returns a dict with
    variant_id, consequence, gene, hgvsp, and per-pop af. Returns None
    only on hard error; returns an empty dict {} (cached) when rsID
    just isn't in gnomAD."""
    if not rsid or not rsid.startswith("rs"):
        return None
    cache = cache_dir / "gnomad" / f"{rsid}.json"
    if use_cache and cache.exists():
        try:
            return json.loads(cache.read_text())
        except Exception:
            pass
    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(
            GNOMAD_URL, method="POST",
            headers={"Content-Type": "application/json", "User-Agent": UA},
            data=json.dumps({"query": _GNOMAD_QUERY,
                              "variables": {"rsid": rsid}}).encode())
        with urllib.request.urlopen(req, timeout=15) as r:
            payload = json.loads(r.read())
    except Exception as e:
        sys.stderr.write(f"[variant-annotate] gnomAD fetch failed for "
                          f"{rsid}: {type(e).__name__}: {e}\n")
        return None
    v = (payload.get("data", {}) or {}).get("variant") or {}
    if not v:
        cache.write_text("{}"); return {}
    # Distill: most-severe consequence + gene + HGVS protein.
    best_csq, best_gene, best_hgvsp, best_rank = None, None, None, 99
    for tc in (v.get("transcript_consequences") or []):
        for term in (tc.get("consequence_terms") or []):
            rk = CSQ_SEVERITY.get(term, 99)
            if rk < best_rank:
                best_rank, best_csq = rk, term
                best_gene = tc.get("gene_symbol")
                best_hgvsp = tc.get("hgvsp") or tc.get("hgvsc")
    # Allele frequencies — prefer exome over genome; af = ac/an.
    af_by_pop: dict = {}
    for src in ("exome", "genome"):
        block = v.get(src) or {}
        for pop in (block.get("populations") or []):
            pid = pop.get("id") or ""
            if pid and "_" not in pid and pid not in af_by_pop:
                ac, an = pop.get("ac"), pop.get("an")
                af_by_pop[pid.lower()] = (ac / an) if (ac is not None
                                                        and an) else None
    out = {
        "variant_id": v.get("variant_id"),
        "consequence": best_csq,
        "gene": best_gene,
        "hgvsp": best_hgvsp,
        "af": af_by_pop,
    }
    cache.write_text(json.dumps(out))
    return out


def gtex_query(rsid: str, cache_dir: Path,
                use_cache: bool = True) -> Optional[dict]:
    """Fetch GTEx v8 top eQTL + sQTL for an rsID. Returns a dict with
    gtex_variant_id, top_eqtl, top_sqtl (each tissue/gene/nes/p), or
    None on hard failure."""
    if not rsid or not rsid.startswith("rs"):
        return None
    cache = cache_dir / "gtex" / f"{rsid}.json"
    if use_cache and cache.exists():
        try:
            return json.loads(cache.read_text())
        except Exception:
            pass
    cache.parent.mkdir(parents=True, exist_ok=True)
    out: dict = {"gtex_variant_id": None, "top_eqtl": None, "top_sqtl": None,
                  "n_eqtl_tissues": 0, "n_sqtl_tissues": 0}
    # Step 1: rsid → GTEx variantId.
    try:
        req = urllib.request.Request(
            f"{GTEX_BASE}/dataset/variant?snpId={rsid}",
            headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=15) as r:
            payload = json.loads(r.read())
        rows = payload.get("data") or []
        if not rows:
            cache.write_text(json.dumps(out)); return out
        out["gtex_variant_id"] = rows[0].get("variantId")
    except Exception as e:
        sys.stderr.write(f"[variant-annotate] GTEx variant lookup failed "
                          f"for {rsid}: {type(e).__name__}: {e}\n")
        return None
    # Step 2: pull single-tissue eQTL + sQTL, pick the smallest-p row.
    gtex_id = out["gtex_variant_id"]
    for endpoint, top_key, n_key in (
            ("singleTissueEqtl", "top_eqtl", "n_eqtl_tissues"),
            ("singleTissueSqtl", "top_sqtl", "n_sqtl_tissues")):
        try:
            url = (f"{GTEX_BASE}/association/{endpoint}"
                   f"?variantId={gtex_id}&itemsPerPage=250")
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=20) as r:
                payload = json.loads(r.read())
            rows = [row for row in (payload.get("data") or [])
                     if row.get("pValue") is not None]
            if not rows:
                continue
            rows.sort(key=lambda row: float(row.get("pValue") or 1))
            top = rows[0]
            out[top_key] = {
                "tissue": top.get("tissueSiteDetailId"),
                "gene": top.get("geneSymbol"),
                "nes": top.get("nes"),
                "p": top.get("pValue"),
            }
            out[n_key] = len({row.get("tissueSiteDetailId") for row in rows})
        except Exception as e:
            sys.stderr.write(f"[variant-annotate] GTEx {endpoint} failed "
                              f"for {rsid}: {type(e).__name__}: {e}\n")
    cache.write_text(json.dumps(out))
    return out


def annotate_rsids(rsids: list[str], sources: set[str],
                    cache_dir: Path,
                    use_cache: bool = True) -> list[dict]:
    """Annotate a list of rsIDs against the requested sources. Returns
    one dict per input rsID, preserving order."""
    results = []
    for rsid in rsids:
        rec: dict = {"rsid": rsid}
        if "gnomad" in sources:
            g = gnomad_query(rsid, cache_dir, use_cache) or {}
            rec["gnomad"] = g
        if "gtex" in sources:
            t = gtex_query(rsid, cache_dir, use_cache) or {}
            rec["gtex"] = t
        results.append(rec)
    return results


# ----------------------------- Output formats -----------------------------

def _flatten_row(rec: dict, sources: set[str]) -> dict:
    """Flatten a nested {gnomad: {...}, gtex: {...}} record into a flat
    column dict for TSV/markdown output."""
    out = {"rsid": rec.get("rsid")}
    if "gnomad" in sources:
        g = rec.get("gnomad") or {}
        out["gnomad_variant_id"] = g.get("variant_id")
        out["consequence"] = g.get("consequence")
        out["gene"] = g.get("gene")
        out["hgvsp"] = g.get("hgvsp")
        af = g.get("af") or {}
        for pop in GNOMAD_POPS:
            v = af.get(pop)
            out[f"af_{pop}"] = (f"{v:.5g}" if isinstance(v, (int, float))
                                else "")
    if "gtex" in sources:
        t = rec.get("gtex") or {}
        out["gtex_variant_id"] = t.get("gtex_variant_id")
        eq = t.get("top_eqtl") or {}
        sq = t.get("top_sqtl") or {}
        out["top_eqtl_tissue"] = eq.get("tissue")
        out["top_eqtl_gene"] = eq.get("gene")
        out["top_eqtl_nes"] = (f"{eq.get('nes'):.4g}"
                                if isinstance(eq.get("nes"), (int, float))
                                else "")
        out["top_eqtl_p"] = (f"{eq.get('p'):.2e}"
                              if isinstance(eq.get("p"), (int, float))
                              else "")
        out["top_sqtl_tissue"] = sq.get("tissue")
        out["top_sqtl_gene"] = sq.get("gene")
        out["top_sqtl_p"] = (f"{sq.get('p'):.2e}"
                              if isinstance(sq.get("p"), (int, float))
                              else "")
        out["n_eqtl_tissues"] = t.get("n_eqtl_tissues", 0)
        out["n_sqtl_tissues"] = t.get("n_sqtl_tissues", 0)
    return {k: ("" if v is None else v) for k, v in out.items()}


def render_tsv(records: list[dict], sources: set[str], fh) -> None:
    rows = [_flatten_row(r, sources) for r in records]
    if not rows:
        return
    cols = list(rows[0].keys())
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def render_markdown(records: list[dict], sources: set[str], fh) -> None:
    rows = [_flatten_row(r, sources) for r in records]
    if not rows:
        return
    cols = list(rows[0].keys())
    fh.write("| " + " | ".join(cols) + " |\n")
    fh.write("|" + "|".join("---" for _ in cols) + "|\n")
    for r in rows:
        fh.write("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |\n")


def render_json(records: list[dict], sources: set[str], fh) -> None:
    json.dump(records, fh, indent=2); fh.write("\n")


# --------------------------------- CLI ---------------------------------

def _collect_rsids(args) -> list[str]:
    rsids: list[str] = []
    if args.rsid:
        rsids.extend(args.rsid)
    if args.rsid_file:
        with open(args.rsid_file) as f:
            rsids.extend(line.strip() for line in f
                          if line.strip() and not line.startswith("#"))
    if args.tsv:
        if not args.rsid_col:
            sys.exit("ERROR: --tsv requires --rsid-col")
        with open(args.tsv) as f:
            header = f.readline().rstrip("\n").split("\t")
            try:
                idx = header.index(args.rsid_col)
            except ValueError:
                sys.exit(f"ERROR: column {args.rsid_col!r} not in TSV "
                         f"header {header}")
            for line in f:
                cells = line.rstrip("\n").split("\t")
                if idx < len(cells) and cells[idx]:
                    rsids.append(cells[idx])
    # De-duplicate while preserving order.
    seen, out = set(), []
    for r in rsids:
        if r and r not in seen:
            seen.add(r); out.append(r)
    return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="annotate.py",
        description="Per-rsID variant annotation against gnomAD + GTEx.")
    ap.add_argument("--rsid", action="append",
                     help="One rsID. Repeat for multiple.")
    ap.add_argument("--rsid-file",
                     help="File with one rsID per line.")
    ap.add_argument("--tsv",
                     help="TSV file; rsIDs are taken from --rsid-col.")
    ap.add_argument("--rsid-col",
                     help="Column name in --tsv that holds rsIDs.")
    ap.add_argument("--sources", default="both",
                     choices=["gnomad", "gtex", "both"],
                     help="Which annotation sources to fetch (default: both)")
    ap.add_argument("--format", default="tsv",
                     choices=["tsv", "markdown", "json"],
                     help="Output format (default: tsv)")
    ap.add_argument("--out", help="Output path (default: stdout)")
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE),
                     help=f"Per-rsID cache dir (default: {DEFAULT_CACHE})")
    ap.add_argument("--no-cache", action="store_true",
                     help="Skip cache reads; refetch every rsID.")
    args = ap.parse_args(argv)
    rsids = _collect_rsids(args)
    if not rsids:
        sys.exit("ERROR: no rsIDs given (use --rsid, --rsid-file, or "
                 "--tsv + --rsid-col)")
    sources: set[str] = ({"gnomad", "gtex"} if args.sources == "both"
                          else {args.sources})
    cache_dir = Path(args.cache_dir).expanduser()
    sys.stderr.write(f"[variant-annotate] {len(rsids)} rsIDs, sources="
                      f"{','.join(sorted(sources))}, cache={cache_dir}\n")
    records = annotate_rsids(rsids, sources, cache_dir,
                              use_cache=(not args.no_cache))
    renderer = {"tsv": render_tsv, "markdown": render_markdown,
                 "json": render_json}[args.format]
    if args.out:
        with open(args.out, "w") as f:
            renderer(records, sources, f)
        sys.stderr.write(f"[variant-annotate] wrote {args.out}\n")
    else:
        renderer(records, sources, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
