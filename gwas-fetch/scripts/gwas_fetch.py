#!/usr/bin/env python3
"""Fetch harmonised GWAS sumstats from the GWAS Catalog by GCST accession.

Two-tier resolution: OKG-first ($OKG_REPO required for that path), GWAS
Catalog REST API fallback. Writes a sidecar .fetch.json with provenance.
This skill does NOT lift coordinates — chain into the liftover skill for that.

Usage:
    gwas_fetch.py --gcst GCST90704615 [--cache-dir ~/.cache/gwas-catalog]
                  [--refresh] [--okg-repo /path/to/okg]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


BUILD_ALIASES = {
    "hg18": "hg18", "grch36": "hg18",
    "hg19": "hg19", "grch37": "hg19", "b37": "hg19", "hg37": "hg19",
    "hg38": "hg38", "grch38": "hg38", "b38": "hg38",
}


def normalize_build(b: Optional[str]) -> Optional[str]:
    if b is None:
        return None
    return BUILD_ALIASES.get(b.lower().strip())


def _gwas_catalog_ftp_dir(accession: str) -> Optional[str]:
    """Return the FTP base for a GCST accession using the bucket layout
    GCST<lo>-GCST<hi>/<accession>/, where the bucket is the 1000-sized
    window containing the numeric accession (lo=floor((n-1)/1000)*1000+1)."""
    m = re.match(r"GCST(\d+)$", accession)
    if not m:
        return None
    n = int(m.group(1))
    lo = ((n - 1) // 1000) * 1000 + 1
    hi = lo + 999
    return ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            f"GCST{lo:06d}-GCST{hi:06d}/{accession}/")


def _probe_harmonised_url(ftp_dir: str, accession: str) -> Optional[str]:
    """List the harmonised/ subdirectory and pick the .h.tsv.gz file."""
    url = ftp_dir.rstrip("/") + "/harmonised/"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            html = r.read().decode("utf-8", errors="replace")
    except Exception:
        return None
    candidates = [
        m.group(1) for m in re.finditer(r'href="([^"]+\.h\.tsv\.gz)"', html)
    ]
    if not candidates:
        return None
    # Prefer files that mention the accession (defensive against indexing).
    preferred = [c for c in candidates if accession in c]
    chosen = (preferred or candidates)[0]
    return url + chosen


def _probe_root_sumstats_url(ftp_dir: str,
                              accession: str) -> Optional[tuple[str, Optional[str]]]:
    """Fallback for deposits that don't ship a harmonised/ subdirectory.

    Lists the GCST root and picks the first `*build*.tsv.gz` (or any
    `*.tsv.gz`) it finds. Returns `(url, parsed_build)` where parsed_build
    is the build extracted from the filename when present (e.g. `_buildGRCh37`)
    or None. Harmonised lookups should be preferred over this — only fall
    back here when the harmonised path is empty."""
    url = ftp_dir.rstrip("/") + "/"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            html = r.read().decode("utf-8", errors="replace")
    except Exception:
        return None
    candidates = [
        m.group(1) for m in re.finditer(r'href="([^"]+\.tsv\.gz)"', html)
    ]
    # Skip metadata/index files; only keep actual data tarballs.
    candidates = [
        c for c in candidates
        if not c.endswith(("-meta.yaml.gz", ".tbi"))
        and ".tsv.gz" in c
    ]
    if not candidates:
        return None
    preferred = [c for c in candidates if accession in c]
    chosen = (preferred or candidates)[0]
    parsed_build = None
    m = re.search(r"_build(GRCh3[78]|hg1[89])", chosen, re.IGNORECASE)
    if m:
        parsed_build = normalize_build(m.group(1))
    return url + chosen, parsed_build


def gwas_catalog_studies_by_pmid(pmid: str) -> list:
    """Return a list of accessionIds for all GWAS Catalog studies under
    the given PubMed ID."""
    api = ("https://www.ebi.ac.uk/gwas/rest/api/studies/search/"
            f"findByPublicationIdPubmedId?pubmedId={pmid}&size=200")
    try:
        with urllib.request.urlopen(api, timeout=30) as r:
            payload = json.loads(r.read())
    except Exception as e:
        raise RuntimeError(f"PubMed ID lookup failed for {pmid}: {e}")
    studies = payload.get("_embedded", {}).get("studies", [])
    return [s["accessionId"] for s in studies if s.get("accessionId")]


def _looks_like_sumstats_url(url: Optional[str]) -> bool:
    """True if `url` points at a sumstats file (FTP tsv/csv/parquet) and
    not at the catalog landing page or REST endpoint."""
    if not url:
        return False
    u = url.lower()
    if u.endswith((".tsv.gz", ".tsv", ".csv.gz", ".csv",
                    ".parquet", ".vcf.gz", ".h.tsv.gz")):
        return True
    if "/harmonised/" in u or "ftp.ebi.ac.uk" in u:
        return True
    return False


def okg_lookup(accession: str, okg_repo: Path) -> Optional[dict]:
    """Try to find the GCST as an alias in the OKG; return attrs dict if hit."""
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
        print(f"warning: OKG_REPO={okg_repo} has no statgen-analysis "
              f"deployment; skipping OKG lookup", file=sys.stderr)
        return None
    env = os.environ.copy()
    env.setdefault("OKG_DSN",
                   "postgres://postgres:okg@localhost:5449/statgen_analysis")
    proc = subprocess.Popen(
        ["uv", "run", "--extra", "mcp", "python",
         "deployments/statgen-analysis/server.py"],
        cwd=str(okg_repo), env=env,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, bufsize=1,
    )
    def send(m): proc.stdin.write(json.dumps(m) + "\n"); proc.stdin.flush()
    def read(): return json.loads(proc.stdout.readline())
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "gwas-fetch", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": "search",
                         "arguments": {"query": accession, "method": "alias",
                                       "limit": 5}}})
        resp = read()
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    sc = resp.get("result", {}).get("structuredContent") or {}
    for hit in sc.get("results", []):
        # MCP search returns the alias hit at the top level and the
        # resolved node nested under `node`. Some older deployments
        # return node fields at the top level; handle both.
        node = hit.get("node") or hit
        subtype = node.get("subtype")
        if subtype not in ("dataset_metadata", "paper"):
            continue
        a = node.get("attrs") or {}
        return {
            "node_id": node.get("node_id"),
            "subtype": subtype,
            "genome_build": normalize_build(a.get("genome_build")),
            "source_url": a.get("source_url"),
            "provider": a.get("provider"),
            "sample_size": a.get("sample_size") or a.get("initial_sample_size"),
            "ancestry": a.get("ancestry_scope") or a.get("ancestry"),
            "summary_statistics_url": a.get("summary_statistics_url"),
            "accession": a.get("accession"),
        }
    return None


def catalog_api_lookup(accession: str) -> dict:
    """Query the GWAS Catalog REST API for study + sumstats metadata.

    Returns a dict with `genome_build`, `harmonised_url`, plus raw metadata.
    Raises RuntimeError on HTTP failure.
    """
    study_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{accession}"
    # REST is best-effort: it gives us first_author / pmid / build for the
    # manifest but isn't required to construct a download URL. The FTP path
    # is derived from the accession bucket alone, and the build can be
    # parsed from the filename when REST fails or omits it.
    study = {}
    try:
        with urllib.request.urlopen(study_url, timeout=30) as r:
            study = json.loads(r.read())
    except Exception as e:
        print(f"warning: GWAS Catalog REST lookup for {accession} failed "
              f"({type(e).__name__}); falling back to FTP-only resolution",
              file=sys.stderr)
    build = None
    for key in ("summaryStatisticsAssembly", "genomeAssembly",
                 "summary_statistics_assembly"):
        if study.get(key):
            build = normalize_build(study[key]); break
    info = study.get("publicationInfo") or {}
    first_author = (info.get("author") or {}).get("fullname")
    pmid = info.get("pubmedId")
    # GWAS Catalog FTP layout: GCST<bucket-lo>-GCST<bucket-hi>/<accession>/harmonised/
    # where bucket = 1000-sized window around the accession number.
    ftp_dir = _gwas_catalog_ftp_dir(accession)
    # Best-effort URL probe via the summary-statistics API.
    ss_url = None
    try:
        ss_api = (f"https://www.ebi.ac.uk/gwas/summary-statistics/api/"
                  f"studies/{accession}")
        with urllib.request.urlopen(ss_api, timeout=30) as r:
            ss = json.loads(r.read())
        for key in ("summary_statistics_url",
                     "harmonised_summary_statistics_url"):
            if isinstance(ss, dict) and ss.get(key):
                ss_url = ss[key]; break
            if isinstance(ss, dict):
                embedded = ss.get("_embedded", {}).get("studies", [])
                if embedded and embedded[0].get(key):
                    ss_url = embedded[0][key]; break
    except Exception:
        pass
    # Probe the FTP harmonised directory to pick the actual `.h.tsv.gz`.
    if ss_url is None and ftp_dir is not None:
        ss_url = _probe_harmonised_url(ftp_dir, accession)
    # If we resolved a harmonised URL but REST didn't give us a build,
    # default to GRCh38 — the GWAS Catalog harmonisation pipeline always
    # emits files on GRCh38.
    if ss_url is not None and "/harmonised/" in ss_url and build is None:
        build = "hg38"
    # Fallback: some deposits (older, non-harmonised) put the file at
    # <bucket>/<GCST>/GCST*_buildXXX.tsv.gz without a `/harmonised/` dir.
    # Prefer harmonised; only use this when harmonised lookup is empty.
    if ss_url is None and ftp_dir is not None:
        probe = _probe_root_sumstats_url(ftp_dir, accession)
        if probe is not None:
            ss_url, parsed_build = probe
            if build is None and parsed_build:
                build = parsed_build
    return {
        "accession": accession,
        "genome_build": build,
        "first_author": first_author,
        "pmid": pmid,
        "ftp_dir": ftp_dir,
        "harmonised_url": ss_url,
        "study_url": study_url,
        "url_inferred": ss_url is not None and ss_url.endswith(
            f"harmonised/{accession}.h.tsv.gz"),
    }


def emit_coverage_stub(out_dir: Path, accession: str, reason: str) -> Path:
    slug = accession.replace(":", "_").replace("/", "_")
    stub_dir = out_dir / "okg-coverage-stubs" / f"add-dataset-{slug}"
    stub_dir.mkdir(parents=True, exist_ok=True)
    p = stub_dir / "proposal.md"
    p.write_text(f"""# OKG coverage-gap for GWAS Catalog `{accession}`

> Generated by the `gwas-fetch` skill: both the OKG and the GWAS Catalog
> REST API failed to resolve a `genome_build` for this study.

## Why
{reason}

## What Changes
- Add a `dataset_metadata` node id `dataset:gwas_catalog:{accession.lower()}`
  with attrs:
    - `provider: GWAS Catalog`
    - `dataset_id: {accession}`
    - `genome_build: GRCh37` or `GRCh38` (verify via the study record)
    - `source_url: <FTP path to harmonised sumstats>`
    - `access_posture: public_metadata`
- Wire up to the appropriate `paper` and `trait` nodes via `references`.

Scaffold in the okg repo:
```
cd "$OKG_REPO" && openspec new change add-dataset-{slug}
```
""")
    return p


def download(url: str, dest: Path, refresh: bool = False) -> None:
    if dest.exists() and not refresh and dest.stat().st_size > 0:
        print(f"using cached: {dest}", file=sys.stderr); return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url} -> {dest}", file=sys.stderr)
    urllib.request.urlretrieve(url, dest)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gcst", type=str,
                   help="GWAS Catalog accession, e.g. GCST90704615")
    group.add_argument("--pubmed-id", type=str,
                   help="PubMed ID; resolves to all GCSTs under that paper "
                        "via the GWAS Catalog REST API and fetches each.")
    p.add_argument("--cache-dir", type=Path,
                   default=Path.home() / ".cache" / "gwas-catalog")
    p.add_argument("--refresh", action="store_true",
                   help="Re-download even if cached")
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None,
                   help="OKG repo for OKG-first metadata lookup (honors $OKG_REPO)")
    args = p.parse_args()

    args.cache_dir.mkdir(parents=True, exist_ok=True)

    # --pubmed-id: resolve to a list of GCSTs and recurse on each.
    if args.pubmed_id:
        try:
            accessions = gwas_catalog_studies_by_pmid(args.pubmed_id)
        except RuntimeError as e:
            sys.exit(f"REFUSED: {e}")
        if not accessions:
            sys.exit(f"REFUSED: no GWAS Catalog studies found for "
                     f"PubMed ID {args.pubmed_id}.")
        print(f"[pubmed-id {args.pubmed_id}] resolved {len(accessions)} GCSTs: "
              f"{', '.join(accessions)}", file=sys.stderr)
        rc = 0
        for acc in accessions:
            args.gcst = acc
            rc |= _fetch_one(args)
        return rc
    return _fetch_one(args)


def _fetch_one(args) -> int:

    okg_hit = None
    if args.okg_repo is not None:
        okg_hit = okg_lookup(args.gcst, args.okg_repo)

    # Always do an API lookup if we don't already have a harmonised
    # sumstats URL — the OKG dataset_metadata.source_url usually points
    # at the catalog study landing page (HTML), not the sumstats file.
    api_hit = None
    need_api = (
        okg_hit is None
        or okg_hit.get("genome_build") is None
        or not _looks_like_sumstats_url(okg_hit.get("summary_statistics_url"))
    )
    if need_api:
        try:
            api_hit = catalog_api_lookup(args.gcst)
        except RuntimeError as e:
            if okg_hit is None:
                emit_coverage_stub(args.cache_dir, args.gcst, str(e))
                sys.exit(f"REFUSED: {e}")

    # Reconcile build + download URL.
    build = (okg_hit and okg_hit.get("genome_build")) \
        or (api_hit and api_hit.get("genome_build"))
    url = (api_hit and api_hit.get("harmonised_url")) \
        or (okg_hit and _looks_like_sumstats_url(
            okg_hit.get("summary_statistics_url"))
            and okg_hit.get("summary_statistics_url"))
    if not build:
        reason = (f"Neither OKG nor GWAS Catalog REST API yielded a "
                  f"`genome_build` for {args.gcst}.")
        emit_coverage_stub(args.cache_dir, args.gcst, reason)
        sys.exit(f"REFUSED: {reason}")
    if not url:
        reason = (f"No download URL could be constructed for {args.gcst}.")
        emit_coverage_stub(args.cache_dir, args.gcst, reason)
        sys.exit(f"REFUSED: {reason}")

    fname = Path(url).name or f"{args.gcst}.tsv.gz"
    cached = args.cache_dir / fname
    download(url, cached, refresh=args.refresh)
    sha = sha256(cached)

    okg_node_ids = {}
    if okg_hit:
        okg_node_ids[okg_hit["subtype"]] = okg_hit["node_id"]

    manifest = {
        "accession": args.gcst,
        "download_url": url,
        "cached_path": str(cached),
        "sha256": sha,
        "genome_build": build,
        "ancestry": (okg_hit and okg_hit.get("ancestry"))
                     or (api_hit and api_hit.get("ancestry")),
        "sample_size": (okg_hit and okg_hit.get("sample_size")),
        "first_author": api_hit and api_hit.get("first_author"),
        "pmid": api_hit and api_hit.get("pmid"),
        "okg_node_ids": okg_node_ids,
        "provenance": {
            "okg_hit": okg_hit,
            "api_hit": api_hit,
            "url_inferred": api_hit and api_hit.get("url_inferred", False),
        },
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = cached.with_suffix(cached.suffix + ".fetch.json")
    mpath.write_text(json.dumps(manifest, indent=2))

    print(f"fetched {args.gcst}: {cached} ({sha[:12]}...; build={build})")
    print(f"manifest -> {mpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
