#!/usr/bin/env python3
"""LDSC wrapper: munge, heritability (h2), genetic correlation (rg).

Handles install of the CBIIT/ldsc Python 3 / Mac-compatible fork on first
use, downloads the canonical EUR LD-score reference, and writes a sidecar
JSON manifest with OKG provenance when $OKG_REPO is set.

Usage:
    ldsc.py munge --in <tsv> --out <prefix> [--N <int>] [...]
    ldsc.py h2    --in <prefix.sumstats.gz> --out <prefix>
    ldsc.py rg    --in1 <a.sumstats.gz> --in2 <b.sumstats.gz> --out <prefix>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


DEFAULT_REPO = "https://github.com/CBIIT/ldsc"
DEFAULT_LD_SCORES_URL = (
    "https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2"
)
CACHE_ROOT = Path.home() / ".cache" / "ldsc"
LDSC_PY_DEPS = ["numpy", "pandas", "scipy", "bitarray"]


# ---------------------------- Install / cache ----------------------------

def ensure_ldsc_repo(repo_url: str, commit: Optional[str],
                     cache_root: Path, refresh: bool) -> Path:
    """Clone (or git-pull) the LDSC repo into the cache and return its path."""
    repo_dir = cache_root / "repo"
    if refresh and repo_dir.exists():
        shutil.rmtree(repo_dir)
    if not repo_dir.exists():
        repo_dir.parent.mkdir(parents=True, exist_ok=True)
        print(f"cloning {repo_url} -> {repo_dir}", file=sys.stderr)
        subprocess.check_call(
            ["git", "clone", "--depth", "1", repo_url, str(repo_dir)],
            stdout=sys.stderr,
        )
    if commit:
        # Need full history for arbitrary commits — fetch then checkout.
        subprocess.check_call(
            ["git", "fetch", "--unshallow"], cwd=repo_dir, stdout=sys.stderr,
            stderr=sys.stderr,
        )
        subprocess.check_call(
            ["git", "checkout", commit], cwd=repo_dir, stdout=sys.stderr,
        )
    return repo_dir


def get_repo_commit(repo_dir: Path) -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo_dir,
                       capture_output=True, text=True)
    return r.stdout.strip()


def ensure_python_deps() -> None:
    """Install LDSC's required Python deps if any are missing."""
    missing = []
    for dep in LDSC_PY_DEPS:
        try:
            __import__(dep)
        except ImportError:
            missing.append(dep)
    if not missing:
        return
    for cmd in [
        ["uv", "pip", "install", *missing],
        [sys.executable, "-m", "pip", "install", "--user", *missing],
    ]:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                print(f"installed LDSC deps: {missing}", file=sys.stderr)
                return
        except FileNotFoundError:
            continue
    sys.exit(f"ERROR: could not install LDSC deps {missing}; "
             f"`pip install {' '.join(missing)}` manually")


def ensure_ld_scores(ld_scores_dir: Optional[Path],
                     cache_root: Path,
                     url: str = DEFAULT_LD_SCORES_URL,
                     refresh: bool = False) -> tuple[Path, str]:
    """Return (ld_scores_dir, sha256). Downloads + extracts on first use."""
    if ld_scores_dir is not None and ld_scores_dir.exists():
        sha = _sha256_of_dir(ld_scores_dir)
        return ld_scores_dir, sha
    default = cache_root / "ld_scores" / "eur_w_ld_chr"
    if default.exists() and not refresh:
        return default, _sha256_of_dir(default)
    cache_root.mkdir(parents=True, exist_ok=True)
    tarball = cache_root / "eur_w_ld_chr.tar.bz2"
    if refresh and tarball.exists():
        tarball.unlink()
    if not tarball.exists():
        print(f"downloading {url} -> {tarball} (~46 MB)", file=sys.stderr)
        urllib.request.urlretrieve(url, tarball)
    sha = hashlib.sha256()
    with open(tarball, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha.update(chunk)
    extract_root = cache_root / "ld_scores"
    extract_root.mkdir(parents=True, exist_ok=True)
    print(f"extracting {tarball} -> {extract_root}", file=sys.stderr)
    with tarfile.open(tarball, "r:bz2") as t:
        t.extractall(extract_root)
    extracted = extract_root / "eur_w_ld_chr"
    if not extracted.exists():
        # The tarball may extract to a different root; pick the first dir.
        candidates = [p for p in extract_root.iterdir() if p.is_dir()]
        if not candidates:
            sys.exit(f"ERROR: extraction produced no directory in {extract_root}")
        extracted = candidates[0]
    return extracted, sha.hexdigest()


def _sha256_of_dir(d: Path) -> str:
    """Cheap content-hash: SHA-256 over sorted filenames + file sizes."""
    h = hashlib.sha256()
    for p in sorted(d.rglob("*")):
        if p.is_file():
            h.update(str(p.relative_to(d)).encode())
            h.update(str(p.stat().st_size).encode())
    return h.hexdigest()


# ---------------------------- OKG provenance ----------------------------

def resolve_okg_node_ids(okg_repo: Optional[Path]) -> dict:
    """Query the OKG for method/software/paper ldsc node ids. Returns {}
    when $OKG_REPO is unset or the OKG MCP can't be reached."""
    if okg_repo is None:
        return {}
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
        print(f"warning: OKG_REPO={okg_repo} doesn't have a statgen-analysis "
              f"server; skipping OKG provenance", file=sys.stderr)
        return {}
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
    out = {}
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "ldsc-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        for nid in ("method:ldsc", "software:ldsc", "paper:ldsc_2015"):
            send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                  "params": {"name": "get_node",
                             "arguments": {"node_id": nid}}})
            resp = read()
            sc = resp.get("result", {}).get("structuredContent") or {}
            if sc.get("node_id") == nid:
                kind = nid.split(":", 1)[0]
                out[kind] = nid
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return out


# ---------------------------- Subcommands ----------------------------

def run_munge(args, repo_dir: Path, okg_node_ids: dict) -> int:
    cmd = [sys.executable, str(repo_dir / "munge_sumstats.py"),
           "--sumstats", str(args.input),
           "--out", str(args.out)]
    if args.N is not None:
        cmd.extend(["--N", str(args.N)])
    if args.N_col is not None:
        cmd.extend(["--N-col", args.N_col])
    if args.snp_col is not None:
        cmd.extend(["--snp", args.snp_col])
    if args.a1_col is not None:
        cmd.extend(["--a1", args.a1_col])
    if args.a2_col is not None:
        cmd.extend(["--a2", args.a2_col])
    if args.p_col is not None:
        cmd.extend(["--p", args.p_col])
    if args.frq_col is not None:
        cmd.extend(["--frq", args.frq_col])
    if args.signed_sumstats:
        cmd.extend(["--signed-sumstats", args.signed_sumstats])
    if args.extra:
        cmd.extend(args.extra)
    print(f"[ldsc munge] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)
    log = Path(str(args.out) + ".log")
    summary = _parse_munge_log(log) if log.exists() else {}
    _write_manifest(args.out, "munge", args, repo_dir, None, summary,
                     okg_node_ids)
    return rc


def run_h2(args, repo_dir: Path, ld_scores_dir: Path, ld_sha: str,
           okg_node_ids: dict) -> int:
    cmd = [sys.executable, str(repo_dir / "ldsc.py"),
           "--h2", str(args.input),
           "--ref-ld-chr", str(ld_scores_dir) + "/",
           "--w-ld-chr", str(ld_scores_dir) + "/",
           "--out", str(args.out)]
    if args.extra:
        cmd.extend(args.extra)
    print(f"[ldsc h2] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)
    log = Path(str(args.out) + ".log")
    summary = _parse_h2_log(log) if log.exists() else {}
    _write_manifest(args.out, "h2", args, repo_dir,
                     {"ld_scores_dir": str(ld_scores_dir),
                      "ld_scores_sha256": ld_sha},
                     summary, okg_node_ids)
    return rc


def run_rg(args, repo_dir: Path, ld_scores_dir: Path, ld_sha: str,
           okg_node_ids: dict) -> int:
    cmd = [sys.executable, str(repo_dir / "ldsc.py"),
           "--rg", f"{args.in1},{args.in2}",
           "--ref-ld-chr", str(ld_scores_dir) + "/",
           "--w-ld-chr", str(ld_scores_dir) + "/",
           "--out", str(args.out)]
    if args.extra:
        cmd.extend(args.extra)
    print(f"[ldsc rg] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)
    log = Path(str(args.out) + ".log")
    summary = _parse_rg_log(log) if log.exists() else {}
    _write_manifest(args.out, "rg", args, repo_dir,
                     {"ld_scores_dir": str(ld_scores_dir),
                      "ld_scores_sha256": ld_sha,
                      "input_1": str(args.in1), "input_2": str(args.in2)},
                     summary, okg_node_ids)
    return rc


# ---------------------------- Log parsers ----------------------------

def _parse_munge_log(log: Path) -> dict:
    text = log.read_text(errors="ignore")
    summary = {}
    for pat, key, cast in [
        (r"(\d+) SNPs remain", "n_snps_remain", int),
        (r"Mean chi\^2 = ([\d.eE+-]+)", "mean_chi2", float),
        (r"Removed (\d+) SNPs with duplicated rs numbers", "n_dup", int),
    ]:
        m = re.search(pat, text)
        if m:
            try: summary[key] = cast(m.group(1))
            except Exception: pass
    return summary


def _parse_h2_log(log: Path) -> dict:
    text = log.read_text(errors="ignore")
    summary = {}
    m = re.search(r"Total Observed scale h2:\s+([-\d.eE+]+)\s+\(([-\d.eE+]+)\)", text)
    if m:
        summary["h2"] = float(m.group(1)); summary["h2_se"] = float(m.group(2))
    m = re.search(r"Intercept:\s+([-\d.eE+]+)\s+\(([-\d.eE+]+)\)", text)
    if m:
        summary["intercept"] = float(m.group(1))
        summary["intercept_se"] = float(m.group(2))
    m = re.search(r"Ratio:\s+([-\d.eE+]+)\s+\(([-\d.eE+]+)\)", text)
    if m:
        summary["ratio"] = float(m.group(1)); summary["ratio_se"] = float(m.group(2))
    m = re.search(r"Mean Chi\^2:\s+([\d.eE+]+)", text)
    if m:
        summary["mean_chi2"] = float(m.group(1))
    return summary


def _parse_rg_log(log: Path) -> dict:
    text = log.read_text(errors="ignore")
    summary = {}
    m = re.search(r"Genetic Correlation:\s+([-\d.eE+]+)\s+\(([-\d.eE+]+)\)", text)
    if m:
        summary["rg"] = float(m.group(1)); summary["rg_se"] = float(m.group(2))
    m = re.search(r"P:\s+([\d.eE+-]+)", text)
    if m:
        summary["p_value"] = float(m.group(1))
    return summary


# ---------------------------- Manifest ----------------------------

def _write_manifest(out_prefix: Path, subcmd: str, args, repo_dir: Path,
                    ld_info: Optional[dict], summary: dict,
                    okg_node_ids: dict) -> Path:
    manifest = {
        "subcommand": subcmd,
        "output_prefix": str(out_prefix),
        "ldsc_repo": str(args.repo_url),
        "ldsc_commit": get_repo_commit(repo_dir),
        "key_results": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    if subcmd == "munge":
        manifest["input"] = str(args.input)
    elif subcmd == "h2":
        manifest["input"] = str(args.input)
    elif subcmd == "rg":
        manifest["input_1"] = str(args.in1)
        manifest["input_2"] = str(args.in2)
    if ld_info:
        manifest.update(ld_info)
    mpath = Path(str(out_prefix) + ".ldsc.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    return mpath


# ---------------------------- CLI ----------------------------

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--repo-url", type=str, default=DEFAULT_REPO,
                   help=f"LDSC repo URL (default: {DEFAULT_REPO})")
    p.add_argument("--repo-commit", type=str, default=None,
                   help="Pin to a specific commit (default: main HEAD)")
    p.add_argument("--repo-cache", type=Path, default=CACHE_ROOT,
                   help=f"Cache root (default: {CACHE_ROOT})")
    p.add_argument("--ld-scores-dir", type=Path, default=None,
                   help="Path to an eur_w_ld_chr/-style dir; default: auto-download")
    p.add_argument("--refresh", action="store_true",
                   help="Force re-clone of LDSC + re-download of LD scores")
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None,
                   help="OKG repo for provenance lookup (honors $OKG_REPO).")
    p.add_argument("extra", nargs=argparse.REMAINDER,
                   help="Extra flags forwarded verbatim to the LDSC subcommand")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="subcmd", required=True)

    pm = sub.add_parser("munge", help="Munge raw sumstats into LDSC format")
    pm.add_argument("--in", dest="input", type=Path, required=True)
    pm.add_argument("--out", type=Path, required=True,
                    help="Output prefix (LDSC writes <out>.sumstats.gz + .log)")
    pm.add_argument("--N", type=int, default=None,
                    help="Constant sample size if not in a column")
    pm.add_argument("--N-col", type=str, default=None)
    pm.add_argument("--snp-col", type=str, default=None)
    pm.add_argument("--a1-col", type=str, default=None,
                    help="Effect allele column name")
    pm.add_argument("--a2-col", type=str, default=None,
                    help="Other allele column name")
    pm.add_argument("--p-col", type=str, default=None)
    pm.add_argument("--frq-col", type=str, default=None,
                    help="Effect allele frequency column")
    pm.add_argument("--signed-sumstats", type=str, default=None,
                    help="LDSC's --signed-sumstats flag value (e.g. 'beta,0')")
    _add_common(pm)

    ph = sub.add_parser("h2", help="Estimate SNP heritability")
    ph.add_argument("--in", dest="input", type=Path, required=True,
                    help="Munged sumstats (.sumstats.gz)")
    ph.add_argument("--out", type=Path, required=True)
    _add_common(ph)

    prg = sub.add_parser("rg", help="Estimate genetic correlation")
    prg.add_argument("--in1", type=Path, required=True)
    prg.add_argument("--in2", type=Path, required=True)
    prg.add_argument("--out", type=Path, required=True)
    _add_common(prg)

    args = p.parse_args()

    # Common setup for every subcommand.
    repo_dir = ensure_ldsc_repo(args.repo_url, args.repo_commit,
                                  Path(args.repo_cache), args.refresh)
    ensure_python_deps()
    okg_node_ids = resolve_okg_node_ids(args.okg_repo)

    if args.subcmd == "munge":
        return run_munge(args, repo_dir, okg_node_ids)
    ld_dir, ld_sha = ensure_ld_scores(args.ld_scores_dir,
                                       Path(args.repo_cache),
                                       refresh=args.refresh)
    if args.subcmd == "h2":
        return run_h2(args, repo_dir, ld_dir, ld_sha, okg_node_ids)
    if args.subcmd == "rg":
        return run_rg(args, repo_dir, ld_dir, ld_sha, okg_node_ids)
    p.error(f"unknown subcommand: {args.subcmd}")


if __name__ == "__main__":
    sys.exit(main())
