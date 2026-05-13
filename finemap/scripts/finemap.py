#!/usr/bin/env python3
"""Fine-mapping wrapper around mancusolab/sushie.

Subcommands:
    finemap.py susie     --vcf <one.vcf> --pheno <one.pheno> [--covar ...] --out <prefix>
    finemap.py sushie    --vcf <k.vcf ...> --pheno <k.pheno ...> [--covar ...] --out <prefix>
    finemap.py sumstats  --z <k.z.tsv ...> --ld <k.ld ...> --n <n1 n2 ...> --out <prefix>
    finemap.py --verify-install   # run the bundled 3-ancestry tutorial

On first use the script clones https://github.com/mancusolab/sushie to
~/.cache/sushie/repo and pip-installs it, then runs the bundled tutorial
under <repo>/data/ as a smoke test. Records OKG provenance to a
.finemap.json sidecar when $OKG_REPO is set.

The `sumstats` mode calls sushie's `infer_sushie_ss` programmatically and
works for K=1 (single-ancestry, equivalent to SuSiE-RSS) or K>=2
(multi-ancestry SuShiE) — the choice is implicit in how many --z / --ld
files you pass.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


CACHE_ROOT = Path.home() / ".cache" / "sushie"
DEFAULT_REPO = "https://github.com/mancusolab/sushie"

OKG_NODES_SUSHIE = {
    "method": "method:sushie",
    "software": "software:sushie",
    "paper": "paper:sushie_2025",
}
OKG_NODES_SUSIE = {
    "method": "method:susie_finemapping",
    "software_operational": "software:sushie",
    "software_canonical": "software:susie",
    "paper_canonical": "paper:susie_2020",
}


# ---------------------------- Install ----------------------------

def ensure_repo(repo_url: str, commit: Optional[str], cache_root: Path,
                refresh: bool = False) -> Path:
    repo_dir = cache_root / "repo"
    if refresh and repo_dir.exists():
        shutil.rmtree(repo_dir)
    if not repo_dir.exists():
        cache_root.mkdir(parents=True, exist_ok=True)
        print(f"cloning {repo_url} -> {repo_dir}", file=sys.stderr)
        subprocess.check_call(
            ["git", "clone", "--depth", "1", repo_url, str(repo_dir)],
            stdout=sys.stderr,
        )
    if commit:
        subprocess.check_call(
            ["git", "fetch", "--unshallow"], cwd=repo_dir,
            stdout=sys.stderr, stderr=sys.stderr,
        )
        subprocess.check_call(
            ["git", "checkout", commit], cwd=repo_dir, stdout=sys.stderr,
        )
    return repo_dir


def ensure_sushie_installed(repo_dir: Path) -> None:
    try:
        import sushie  # noqa: F401
        return
    except ImportError:
        pass
    for cmd in [
        ["uv", "pip", "install", str(repo_dir)],
        [sys.executable, "-m", "pip", "install", "--user", str(repo_dir)],
        [sys.executable, "-m", "pip", "install", str(repo_dir)],
    ]:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                print(f"installed sushie via: {' '.join(cmd[:3])}",
                      file=sys.stderr)
                return
        except FileNotFoundError:
            continue
    sys.exit("ERROR: could not pip-install sushie from "
             f"{repo_dir}; try `pip install {repo_dir}` manually")


def get_repo_commit(repo_dir: Path) -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo_dir,
                       capture_output=True, text=True)
    return r.stdout.strip()


def which_sushie_cli() -> str:
    """Return the path to the sushie CLI, prefer `sushie` on PATH else
    `python -m sushie` as fallback."""
    p = shutil.which("sushie")
    return p if p else f"{sys.executable} -m sushie"


# ---------------------------- Verify install ----------------------------

def verify_install(repo_dir: Path) -> int:
    """Run the bundled 3-ancestry tutorial to verify the install works."""
    sentinel = CACHE_ROOT / "verify_install" / ".verified"
    data_dir = repo_dir / "data"
    out_dir = CACHE_ROOT / "verify_install"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "test_result"
    if not data_dir.exists():
        sys.exit(f"ERROR: bundled tutorial data not found at {data_dir}")
    cli = which_sushie_cli()
    cmd = (cli.split() + [
        "finemap",
        "--pheno", "EUR.pheno", "AFR.pheno", "EAS.pheno",
        "--vcf",
        str(data_dir / "vcf" / "EUR.vcf"),
        str(data_dir / "vcf" / "AFR.vcf"),
        str(data_dir / "vcf" / "EAS.vcf"),
        "--covar", "EUR.covar", "AFR.covar", "EAS.covar",
        "--output", str(out_prefix),
    ])
    print(f"[verify-install] running 3-ancestry tutorial in {data_dir}",
          file=sys.stderr)
    print(f"  cmd: {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd, cwd=str(data_dir))
    if rc != 0:
        sys.exit(f"ERROR: sushie finemap tutorial failed (exit {rc}); "
                 f"see logs above. Sentinel NOT written.")
    # Check that the expected output exists.
    found_any = any(out_dir.glob("test_result*"))
    if not found_any:
        sys.exit(f"ERROR: sushie finemap returned 0 but no test_result* "
                 f"files were created under {out_dir}")
    sentinel.write_text(json.dumps({
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "tutorial": "3-ancestry EUR+AFR+EAS HapMap",
        "outputs": [str(p) for p in sorted(out_dir.glob("test_result*"))],
    }, indent=2))
    print(f"[verify-install] OK; sentinel -> {sentinel}", file=sys.stderr)
    return 0


def is_verified() -> bool:
    return (CACHE_ROOT / "verify_install" / ".verified").exists()


# ---------------------------- OKG resolver ----------------------------

def resolve_okg(okg_repo: Optional[Path], node_ids: list[str]) -> dict:
    if okg_repo is None:
        return {}
    if not (okg_repo / "deployments/statgen-analysis/server.py").exists():
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
    found = {}
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "finemap-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        for nid in node_ids:
            send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                  "params": {"name": "get_node",
                             "arguments": {"node_id": nid}}})
            resp = read()
            sc = resp.get("result", {}).get("structuredContent") or {}
            if sc.get("node_id") == nid:
                found[nid] = nid
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return found


# ---------------------------- Run sushie ----------------------------

def run_finemap(args, mode: str, repo_dir: Path) -> int:
    """Invoke `sushie finemap` with the user's flags. Mode is 'susie' or
    'sushie' — chosen by number of VCF args."""
    cli = which_sushie_cli().split()
    cmd = cli + ["finemap",
                 "--vcf", *[str(p) for p in args.vcf],
                 "--pheno", *[str(p) for p in args.pheno],
                 "--output", str(args.out)]
    if args.covar:
        cmd += ["--covar", *[str(p) for p in args.covar]]
    if args.extra:
        cmd += args.extra
    print(f"[finemap {mode}] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)

    # OKG provenance.
    if mode == "sushie":
        nids = list(OKG_NODES_SUSHIE.values())
        found = resolve_okg(args.okg_repo, nids)
        okg_node_ids = {k: v for k, v in OKG_NODES_SUSHIE.items()
                         if v in found}
    else:
        nids = list(OKG_NODES_SUSIE.values())
        found = resolve_okg(args.okg_repo, nids)
        okg_node_ids = {k: v for k, v in OKG_NODES_SUSIE.items()
                         if v in found}

    # Summary from any test_result.cs.tsv style output.
    summary = _parse_outputs(args.out)
    manifest = {
        "subcommand": mode,
        "output_prefix": str(args.out),
        "sushie_repo": args.repo_url,
        "sushie_commit": get_repo_commit(repo_dir),
        "n_ancestries": len(args.vcf),
        "inputs": {
            "vcf":   [str(p) for p in args.vcf],
            "pheno": [str(p) for p in args.pheno],
            "covar": [str(p) for p in (args.covar or [])],
        },
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(str(args.out) + ".finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    if rc == 0 and summary:
        print(f"summary: {summary.get('n_cs', '?')} CSs, "
              f"mean size {summary.get('mean_cs_size', '?')}, "
              f"max PIP {summary.get('max_pip', '?')}")
    return rc


def _parse_outputs(out_prefix: Path) -> dict:
    summary: dict = {}
    cs_path = Path(f"{out_prefix}.cs.tsv")
    if cs_path.exists():
        try:
            with open(cs_path) as f:
                header = f.readline().rstrip("\n").split("\t")
                rows = [line.rstrip("\n").split("\t") for line in f]
            if rows and "CSIndex" in header:
                idx = header.index("CSIndex")
                cs_groups = {}
                for r in rows:
                    cs_groups.setdefault(r[idx], []).append(r)
                summary["n_cs"] = len(cs_groups)
                if cs_groups:
                    sizes = [len(v) for v in cs_groups.values()]
                    summary["mean_cs_size"] = round(sum(sizes) / len(sizes), 2)
                    summary["cs_sizes"] = sorted(sizes)
        except Exception as e:
            summary["cs_parse_error"] = str(e)
    weight_path = Path(f"{out_prefix}.weight.tsv")
    if weight_path.exists():
        try:
            with open(weight_path) as f:
                header = f.readline().rstrip("\n").split("\t")
                pip_idx = header.index("PIP") if "PIP" in header else None
                if pip_idx is not None:
                    max_pip = 0.0
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        try:
                            v = float(parts[pip_idx])
                            if v > max_pip:
                                max_pip = v
                        except (ValueError, IndexError):
                            continue
                    summary["max_pip"] = round(max_pip, 4)
        except Exception as e:
            summary["weight_parse_error"] = str(e)
    corr_path = Path(f"{out_prefix}.corr.tsv")
    if corr_path.exists():
        summary["cross_ancestry_correlation_file"] = str(corr_path)
    return summary


# ---------------------------- Sumstats mode ----------------------------

def _load_z(path: Path) -> tuple[list[str], list[float]]:
    """Load a z-score file. Accept TSV/CSV/whitespace with a header. If a
    column named Z/z is present, use it. Else compute z = beta / se from
    BETA/beta + SE/se. Returns (snp_ids, z) — snp_ids may be index strings
    if no SNP/snp/rsid/SNPID column is found."""
    import csv
    with open(path) as f:
        first = f.readline()
    sep = "\t" if "\t" in first else (
        "," if "," in first and ";" not in first else None)
    rows = []
    with open(path) as f:
        if sep:
            reader = csv.reader(f, delimiter=sep)
        else:
            reader = ([x for x in line.split()] for line in f)
        header = next(iter(reader))
        for r in reader:
            if r and any(x.strip() for x in r):
                rows.append(r)
    cols = {c.lower(): i for i, c in enumerate(header)}
    snp_col = next((cols[k] for k in ("snp", "snpid", "rsid", "snp_id",
                                       "variant_id", "id") if k in cols), None)
    z_col = next((cols[k] for k in ("z", "zscore", "z_score") if k in cols),
                  None)
    if z_col is not None:
        snps = [r[snp_col] if snp_col is not None else str(i)
                for i, r in enumerate(rows)]
        zs = [float(r[z_col]) for r in rows]
        return snps, zs
    beta_col = next((cols[k] for k in ("beta", "effect_size", "b")
                      if k in cols), None)
    se_col = next((cols[k] for k in ("se", "standard_error", "stderr")
                    if k in cols), None)
    if beta_col is None or se_col is None:
        sys.exit(f"ERROR: {path} has no Z/z column and no BETA+SE columns "
                 f"to compute it; header was {header}")
    snps = [r[snp_col] if snp_col is not None else str(i)
            for i, r in enumerate(rows)]
    zs = [float(r[beta_col]) / float(r[se_col]) for r in rows]
    return snps, zs


def _load_ld(path: Path):
    """Load an LD matrix. Supports .npy (numpy save) and whitespace TSV."""
    import numpy as np
    p = Path(path)
    if p.suffix == ".npy":
        ld = np.load(p)
    else:
        ld = np.loadtxt(p)
    if ld.ndim != 2 or ld.shape[0] != ld.shape[1]:
        sys.exit(f"ERROR: LD at {p} is not a square matrix "
                 f"(shape={ld.shape})")
    return ld


def run_sumstats(args, repo_dir: Path) -> int:
    """Programmatic SuShiE/SuSiE sumstats fine-mapping via infer_sushie_ss."""
    import numpy as np
    try:
        from sushie.infer_ss import infer_sushie_ss
    except ImportError:
        sys.exit("ERROR: sushie package not importable; "
                 "run --verify-install first")

    K = len(args.z)
    if len(args.ld) != K:
        sys.exit(f"ERROR: --z count ({K}) != --ld count ({len(args.ld)})")
    if len(args.n) != K:
        sys.exit(f"ERROR: --z count ({K}) != --n count ({len(args.n)})")

    # Load + align by SNP if all files carry SNP IDs; else by row order.
    z_arrs, snp_lists, lds = [], [], []
    for i in range(K):
        snps_i, zs_i = _load_z(Path(args.z[i]))
        ld_i = _load_ld(Path(args.ld[i]))
        if ld_i.shape[0] != len(zs_i):
            sys.exit(f"ERROR: ancestry {i+1}: z file has {len(zs_i)} rows "
                     f"but LD is {ld_i.shape[0]}x{ld_i.shape[0]}")
        z_arrs.append(np.asarray(zs_i, dtype=np.float64))
        snp_lists.append(snps_i)
        lds.append(ld_i.astype(np.float64))

    # Sanity check: all ancestries must have same number of SNPs.
    ms = {len(z) for z in z_arrs}
    if len(ms) != 1:
        sys.exit(f"ERROR: ancestries have mismatched SNP counts: "
                 f"{[len(z) for z in z_arrs]}")
    m = next(iter(ms))

    # Run infer_sushie_ss.
    ns = np.asarray([int(x) for x in args.n], dtype=np.float64)
    print(f"[sumstats] K={K} ancestries, m={m} SNPs, L={args.L}, "
          f"max_iter={args.max_iter}", file=sys.stderr)
    result = infer_sushie_ss(
        lds=np.asarray(lds), ns=ns, zs=np.asarray(z_arrs),
        L=args.L, max_iter=args.max_iter, min_tol=args.min_tol,
        threshold=args.threshold, purity=args.purity,
        min_snps=max(10, args.min_snps),
    )

    # Write outputs in a finemap-friendly format.
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cs_path = Path(f"{out_prefix}.cs.tsv")
    weight_path = Path(f"{out_prefix}.weight.tsv")

    cs_df = result.cs  # pd.DataFrame with SNPIndex, CSIndex, alpha, c_alpha
    snps0 = snp_lists[0]
    if hasattr(cs_df, "to_csv"):
        # Resolve SNP IDs from ancestry 1 (all should match).
        cs_out = cs_df.copy()
        try:
            cs_out["SNP"] = [snps0[int(i)] for i in cs_out["SNPIndex"].values]
        except Exception:
            pass
        cs_out.to_csv(cs_path, sep="\t", index=False)
    pip_arr = np.asarray(result.pip)
    with open(weight_path, "w") as f:
        f.write("SNPIndex\tSNP\tPIP\n")
        for i in range(m):
            f.write(f"{i}\t{snps0[i]}\t{pip_arr[i]:.6g}\n")

    # OKG provenance.
    if K == 1:
        node_map = OKG_NODES_SUSIE
    else:
        node_map = OKG_NODES_SUSHIE
    found = resolve_okg(args.okg_repo, list(node_map.values()))
    okg_node_ids = {k: v for k, v in node_map.items() if v in found}

    summary = _parse_outputs(out_prefix)
    if "max_pip" not in summary:
        summary["max_pip"] = round(float(pip_arr.max()), 4)
    manifest = {
        "subcommand": "sumstats",
        "output_prefix": str(out_prefix),
        "sushie_repo": args.repo_url,
        "sushie_commit": get_repo_commit(repo_dir),
        "n_ancestries": K,
        "n_snps": m,
        "inputs": {
            "z":  [str(p) for p in args.z],
            "ld": [str(p) for p in args.ld],
            "n":  [int(x) for x in args.n],
        },
        "infer_params": {
            "L": args.L, "max_iter": args.max_iter,
            "min_tol": args.min_tol, "threshold": args.threshold,
            "purity": args.purity, "min_snps": args.min_snps,
        },
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{out_prefix}.finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    print(f"summary: {summary.get('n_cs', '?')} CSs, max PIP "
          f"{summary.get('max_pip', '?')}")
    return 0


# ---------------------------- CLI ----------------------------

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--vcf", type=Path, nargs="+", required=True,
                   help="One or more VCF files (one per ancestry)")
    p.add_argument("--pheno", type=Path, nargs="+", required=True,
                   help="One or more phenotype files (sushie format)")
    p.add_argument("--covar", type=Path, nargs="+", default=None)
    p.add_argument("--out", type=Path, required=True,
                   help="Output prefix")
    p.add_argument("--repo-url", type=str, default=DEFAULT_REPO)
    p.add_argument("--sushie-commit", dest="commit", type=str, default=None,
                   help="Pin sushie to a specific commit")
    p.add_argument("--repo-cache", type=Path, default=CACHE_ROOT)
    p.add_argument("--refresh", action="store_true",
                   help="Re-clone sushie and re-run verify-install")
    p.add_argument("--re-verify", action="store_true",
                   help="Force re-running verify-install even if sentinel exists")
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None,
                   help="OKG repo path (honors $OKG_REPO)")
    p.add_argument("extra", nargs=argparse.REMAINDER,
                   help="Extra flags forwarded to `sushie finemap`")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--verify-install", action="store_true",
                   help="Run the bundled 3-ancestry tutorial and exit")
    sub = p.add_subparsers(dest="subcmd")

    psusie = sub.add_parser("susie", help="Single-ancestry fine-mapping")
    _add_common(psusie)
    psushie = sub.add_parser("sushie", help="Multi-ancestry fine-mapping")
    _add_common(psushie)

    pss = sub.add_parser("sumstats",
                          help="Sumstats fine-mapping via infer_sushie_ss "
                               "(K=1 SuSiE-RSS or K>=2 SuShiE)")
    pss.add_argument("--z", type=Path, nargs="+", required=True,
                     help="One or more z-score files (TSV with Z, or "
                          "BETA+SE columns) — one per ancestry")
    pss.add_argument("--ld", type=Path, nargs="+", required=True,
                     help="One or more LD matrices (.npy or whitespace TSV) "
                          "— one per ancestry, same SNP order as --z")
    pss.add_argument("--n", type=int, nargs="+", required=True,
                     help="Sample size per ancestry (one int per --z file)")
    pss.add_argument("--out", type=Path, required=True,
                     help="Output prefix")
    pss.add_argument("--L", type=int, default=10,
                     help="Max number of single effects (default 10)")
    pss.add_argument("--max-iter", type=int, default=500)
    pss.add_argument("--min-tol", type=float, default=1e-4)
    pss.add_argument("--threshold", type=float, default=0.95,
                     help="Credible set coverage threshold (default 0.95)")
    pss.add_argument("--purity", type=float, default=0.5,
                     help="Min CS purity (default 0.5)")
    pss.add_argument("--min-snps", type=int, default=100,
                     help="Min SNPs to fine-map (default 100)")
    pss.add_argument("--repo-url", type=str, default=DEFAULT_REPO)
    pss.add_argument("--sushie-commit", dest="commit", type=str, default=None)
    pss.add_argument("--repo-cache", type=Path, default=CACHE_ROOT)
    pss.add_argument("--refresh", action="store_true")
    pss.add_argument("--re-verify", action="store_true")
    pss.add_argument("--okg-repo", type=Path,
                     default=Path(os.environ["OKG_REPO"])
                             if os.environ.get("OKG_REPO") else None)

    args, _ = p.parse_known_args()

    # Pure verify-install mode.
    if args.verify_install:
        repo_dir = ensure_repo(DEFAULT_REPO, None,
                                CACHE_ROOT, refresh=False)
        ensure_sushie_installed(repo_dir)
        return verify_install(repo_dir)

    if args.subcmd not in ("susie", "sushie", "sumstats"):
        p.print_help(sys.stderr)
        return 2

    # Re-parse for the subcommand's flags.
    args = p.parse_args()

    # Common setup.
    repo_dir = ensure_repo(args.repo_url, args.commit,
                            Path(args.repo_cache), refresh=args.refresh)
    ensure_sushie_installed(repo_dir)

    # First-run verify-install gate.
    if args.refresh or args.re_verify or not is_verified():
        rc = verify_install(repo_dir)
        if rc != 0:
            sys.exit("verify-install failed; refusing to proceed with user run")

    if args.subcmd == "sumstats":
        return run_sumstats(args, repo_dir)

    # Validate ancestry count vs subcommand.
    if args.subcmd == "susie" and len(args.vcf) != 1:
        p.error("--vcf must be a single file in `susie` mode "
                "(got {} files)".format(len(args.vcf)))
    if args.subcmd == "sushie" and len(args.vcf) < 2:
        p.error("`sushie` mode requires K>=2 ancestries; "
                "got {} VCF".format(len(args.vcf)))
    if len(args.vcf) != len(args.pheno):
        p.error("number of --vcf files must match number of --pheno files "
                "({} vs {})".format(len(args.vcf), len(args.pheno)))
    if args.covar and len(args.covar) != len(args.vcf):
        p.error("number of --covar files must match number of --vcf files")

    return run_finemap(args, args.subcmd, repo_dir)


if __name__ == "__main__":
    sys.exit(main())
