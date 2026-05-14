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
    _patch_ldsc_repo(repo_dir)
    return repo_dir


def _patch_ldsc_repo(repo_dir: Path) -> None:
    """Apply small Py3 / pandas-compat patches to the cloned LDSC fork.

    Both upstream `bulik/ldsc` and the CBIIT `Python 3` fork still ship a
    `read_header` that decodes gzipped bytes as if they were `str`. We
    monkey-patch the clone in place so munge_sumstats.py works on .gz
    inputs. Idempotent: applies only if the buggy line is still present.
    """
    ms = repo_dir / "munge_sumstats.py"
    if not ms.exists():
        return
    src = ms.read_text()
    bad = "return [x.rstrip('\\n') for x in openfunc(fh).readline().split()]"
    fix = ("line = openfunc(fh).readline()\n    "
            "if isinstance(line, bytes):\n        line = line.decode('utf-8')\n    "
            "return [x.rstrip('\\n') for x in line.split()]")
    if bad in src and fix.split("\n")[0] not in src:
        ms.write_text(src.replace(bad, fix))
        print(f"  patched read_header in {ms}", file=sys.stderr)


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
                     refresh: bool = False,
                     okg_panel: Optional[dict] = None) -> tuple[Path, str]:
    """Return (ld_scores_dir, sha256). Resolution order:
       1) explicit --ld-scores-dir if it exists,
       2) OKG ld_panel.local_path_hint if it exists locally,
       3) default cache path (~/.cache/ldsc/ld_scores/eur_w_ld_chr),
       4) download from `url` and extract.
    """
    if ld_scores_dir is not None and ld_scores_dir.exists():
        sha = _sha256_of_dir(ld_scores_dir)
        return ld_scores_dir, sha
    if okg_panel and okg_panel.get("local_path_hint"):
        hint = Path(os.path.expanduser(okg_panel["local_path_hint"]))
        if hint.exists():
            print(f"using OKG-resolved LD path: {hint} "
                  f"(from {okg_panel.get('node_id')})", file=sys.stderr)
            return hint, _sha256_of_dir(hint)
    default = cache_root / "ld_scores" / "eur_w_ld_chr"
    if default.exists() and not refresh:
        return default, _sha256_of_dir(default)
    cache_root.mkdir(parents=True, exist_ok=True)
    tarball = cache_root / "eur_w_ld_chr.tar.bz2"
    if refresh and tarball.exists():
        tarball.unlink()
    if not tarball.exists():
        print(f"downloading {url} -> {tarball} (~46 MB)", file=sys.stderr)
        try:
            urllib.request.urlretrieve(url, tarball)
        except Exception as e:
            msg = (f"ERROR: failed to download LD scores from {url}: {e}.")
            if okg_panel and okg_panel.get("source_url"):
                msg += (f"\n  OKG ld_panel {okg_panel.get('node_id')} also "
                        f"points at {okg_panel['source_url']} — that link may "
                        f"require manual download; place the extracted "
                        f"eur_w_ld_chr/ dir at "
                        f"{default} or pass --ld-scores-dir <path>.")
            sys.exit(msg)
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
        # method:ldsc, paper:ldsc_2015: canonical.
        # software:ldsc: upstream bulik/ldsc.
        # software:ldsc_cbiit: the Python-3 / macOS-arm64 fork that this
        #   skill actually installs and runs. Recorded as software_operational
        #   on the manifest so the provenance trail names the running fork.
        for nid in ("method:ldsc", "software:ldsc", "software:ldsc_cbiit",
                     "paper:ldsc_2015"):
            send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                  "params": {"name": "get_node",
                             "arguments": {"node_id": nid}}})
            resp = read()
            sc = resp.get("result", {}).get("structuredContent") or {}
            if sc.get("node_id") == nid:
                if nid == "software:ldsc_cbiit":
                    out["software_operational"] = nid
                elif nid == "software:ldsc":
                    out["software_upstream"] = nid
                else:
                    kind = nid.split(":", 1)[0]
                    out[kind] = nid
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return out


def _strip_dashdash(extras: list) -> list:
    """Drop a literal '--' separator from REMAINDER-captured extras so it
    isn't forwarded into the underlying LDSC CLI (which doesn't expect it)."""
    return [a for a in extras if a != "--"]


def resolve_okg_ld_panel(okg_repo: Optional[Path],
                          panel_id: str) -> Optional[dict]:
    """Query the OKG for an ld_panel node and return its attrs
    (specifically local_path_hint + source_url). Returns None if the OKG
    isn't reachable or the node doesn't exist."""
    if okg_repo is None:
        return None
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
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
                         "clientInfo": {"name": "ldsc-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": "get_node",
                         "arguments": {"node_id": panel_id}}})
        resp = read()
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    sc = resp.get("result", {}).get("structuredContent") or {}
    if sc.get("node_id") != panel_id:
        return None
    a = sc.get("attrs") or {}
    return {
        "node_id": panel_id,
        "local_path_hint": a.get("local_path_hint"),
        "source_url": a.get("source_url"),
        "name": a.get("name"),
        "genome_build": a.get("genome_build"),
    }


# ---------------------------- Sumstats pre-flight ----------------------------

# Columns LDSC's auto-matcher recognises and that the GWAS Catalog harmonised
# files duplicate with a `hm_` prefix. Used by the auto-ignore precheck.
_LDSC_RELEVANT_COLS = (
    "effect_allele", "other_allele", "beta", "odds_ratio",
    "effect_allele_frequency", "p_value", "se", "standard_error",
    "z", "info",
)

_NA_TOKENS = {"NA", "", "nan", "NaN", "N/A", ".", "null", "None"}


def _read_header(path: Path) -> list:
    """Read the header row of a TSV (gz-aware) as a list of column names."""
    import csv
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return next(csv.reader(f, delimiter="\t"))


def _column_is_all_na(path: Path, col: str, sample: int = 50_000) -> bool:
    """Return True if `col` is 100% NA across the first `sample` rows.
    Returns False if the column is missing (let downstream report the
    error) or has any non-NA value."""
    import csv
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return False
        if col not in header:
            return False
        idx = header.index(col)
        for i, row in enumerate(reader):
            if i >= sample:
                break
            if idx >= len(row):
                continue
            if row[idx].strip() not in _NA_TOKENS:
                return False
    return True


def _auto_ignore_non_harmonised(path: Path, chosen_cols: list) -> list:
    """For each user-chosen `hm_*` column whose non-hm sibling exists in the
    same file, return the non-hm name so it can be added to --ignore.
    Also catches LDSC-relevant cols (beta, p_value, etc.) that are present
    in both hm_ and bare forms even when the user didn't explicitly pick the
    hm_ version. Returns deduped list, sorted for determinism."""
    chosen_hm = [c for c in chosen_cols if c and c.startswith("hm_")]
    if not chosen_hm:
        # Without an hm_* signal from the user, don't second-guess: the file
        # might not be a GWAS-Catalog harmonised file at all.
        return []
    try:
        header = _read_header(path)
    except Exception:
        return []
    header_set = set(header)
    ignore = set()
    # User picked hm_X — drop its non-hm sibling if present.
    for hm_col in chosen_hm:
        non_hm = hm_col[3:]
        if non_hm in header_set:
            ignore.add(non_hm)
    # Also catch other LDSC-recognised cols that appear in both forms.
    for col in _LDSC_RELEVANT_COLS:
        if f"hm_{col}" in header_set and col in header_set:
            ignore.add(col)
    return sorted(ignore)


def _merge_ignore_into_extras(extras: list, auto_ignore: list) -> list:
    """Merge auto-detected ignore columns with any --ignore the user passed
    through extras. Returns the new extras list (with --ignore <merged>)."""
    if not auto_ignore:
        return extras
    out = list(extras)
    # Find an existing --ignore and merge.
    for i, tok in enumerate(out):
        if tok == "--ignore" and i + 1 < len(out):
            user_set = {c.strip() for c in out[i + 1].split(",") if c.strip()}
            merged = sorted(user_set | set(auto_ignore))
            out[i + 1] = ",".join(merged)
            return out
    # No existing --ignore — append one.
    out.extend(["--ignore", ",".join(auto_ignore)])
    return out


# ---------------------------- Subcommands ----------------------------

def run_munge(args, repo_dir: Path, okg_node_ids: dict) -> int:
    # No-clobber: if <out>.sumstats.gz exists already, skip the re-munge
    # (idempotent across invocations). Honor --refresh to force.
    out_sumstats = Path(f"{args.out}.sumstats.gz")
    if out_sumstats.exists() and not getattr(args, "refresh", False):
        print(f"[ldsc munge] {out_sumstats} already exists; skipping. "
              f"Pass --refresh to re-munge.", file=sys.stderr)
        log = Path(str(args.out) + ".log")
        summary = _parse_munge_log(log) if log.exists() else {}
        _write_manifest(args.out, "munge", args, repo_dir, None, summary,
                         okg_node_ids)
        return 0

    # Pre-flight: drop --frq-col if the column is 100% NA in this file.
    frq_col = args.frq_col
    if frq_col and getattr(args, "frq_precheck", True):
        if _column_is_all_na(Path(args.input), frq_col):
            print(f"[ldsc munge] warning: column {frq_col!r} is 100% NA in "
                  f"the first 50k rows of {args.input}; dropping --frq-col "
                  f"(LDSC drops every row when the frequency column is NA). "
                  f"Pass --no-frq-precheck to keep it.", file=sys.stderr)
            frq_col = None

    # Pre-flight: auto-add --ignore for non-harmonised duplicates when the
    # user picked hm_* columns from a GWAS-Catalog harmonised file.
    extras = _strip_dashdash(args.extra or [])
    if getattr(args, "auto_ignore_harmonised", True):
        chosen_cols = [args.snp_col, args.a1_col, args.a2_col, args.p_col,
                        frq_col, args.N_col]
        if args.signed_sumstats:
            chosen_cols.append(args.signed_sumstats.split(",", 1)[0])
        auto_ignore = _auto_ignore_non_harmonised(Path(args.input),
                                                    chosen_cols)
        if auto_ignore:
            print(f"[ldsc munge] auto-detected GWAS-Catalog harmonised file; "
                  f"adding --ignore {','.join(auto_ignore)} (non-hm duplicates "
                  f"of LDSC-recognised columns). Pass --no-auto-ignore to "
                  f"disable.", file=sys.stderr)
            extras = _merge_ignore_into_extras(extras, auto_ignore)

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
    if frq_col is not None:
        cmd.extend(["--frq", frq_col])
    if args.signed_sumstats:
        cmd.extend(["--signed-sumstats", args.signed_sumstats])
    if extras:
        cmd.extend(extras)
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
        cmd.extend(_strip_dashdash(args.extra))
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
        cmd.extend(_strip_dashdash(args.extra))
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
    p.add_argument("--okg-ld-panel-id", type=str,
                   default="ld_panel:ldsc_eur_w_ld_chr",
                   help="OKG ld_panel node ID to resolve local_path_hint + "
                        "source_url from (when $OKG_REPO is set). Pass empty "
                        "string to disable.")
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
    pm.add_argument("--no-frq-precheck", dest="frq_precheck",
                    action="store_false",
                    help="Skip the 'is --frq-col all NA?' precheck. By default "
                         "the skill drops --frq-col when the column has no "
                         "non-NA values in the first 50k rows.")
    pm.add_argument("--no-auto-ignore", dest="auto_ignore_harmonised",
                    action="store_false",
                    help="Skip the GWAS-Catalog-harmonised auto --ignore. By "
                         "default the skill adds non-hm duplicates of "
                         "LDSC-recognised columns to --ignore when the user "
                         "chose at least one hm_* column.")
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
    okg_panel = (resolve_okg_ld_panel(args.okg_repo, args.okg_ld_panel_id)
                  if args.okg_ld_panel_id else None)
    ld_dir, ld_sha = ensure_ld_scores(args.ld_scores_dir,
                                       Path(args.repo_cache),
                                       refresh=args.refresh,
                                       okg_panel=okg_panel)
    if okg_panel:
        okg_node_ids["ld_panel"] = okg_panel["node_id"]
    if args.subcmd == "h2":
        return run_h2(args, repo_dir, ld_dir, ld_sha, okg_node_ids)
    if args.subcmd == "rg":
        return run_rg(args, repo_dir, ld_dir, ld_sha, okg_node_ids)
    p.error(f"unknown subcommand: {args.subcmd}")


if __name__ == "__main__":
    sys.exit(main())
