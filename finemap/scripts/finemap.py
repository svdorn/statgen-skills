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
import urllib.request
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


def ensure_sushie_installed(repo_dir: Path) -> Path:
    """Install sushie into an isolated `uv` venv when the host Python is
    incompatible with sushie's pinned deps. Returns the python binary to
    invoke (which `--vcf`/`--gwas`/etc. will then be passed to via
    `python -m sushie`). Reuses the venv on subsequent runs.

    sushie pins glimix_core==3.1.13 which is Python 3.7-3.10 only. On
    Python 3.11+ host we can't install into the host env directly, so
    we use uv to spin up a Python 3.10 venv keyed to the sushie cache.
    """
    venv_dir = repo_dir.parent / "venv-py310"
    venv_py = venv_dir / "bin" / "python"
    if venv_py.exists():
        # Confirm sushie loads.
        r = subprocess.run([str(venv_py), "-c", "import sushie"],
                            capture_output=True, text=True)
        if r.returncode == 0:
            return venv_py

    if not shutil.which("uv"):
        sys.exit("ERROR: `uv` is required to install sushie into a "
                 "Python 3.10 venv (sushie's deps don't support 3.11+). "
                 "Install with `brew install uv` or "
                 "`pipx install uv`.")
    print(f"[sushie] creating Python 3.10 venv at {venv_dir} ...",
          file=sys.stderr)
    r = subprocess.run(["uv", "venv", "--python", "3.10", str(venv_dir)],
                        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: uv venv creation failed: {r.stderr}")
    print(f"[sushie] installing from {repo_dir} into {venv_dir} ...",
          file=sys.stderr)
    r = subprocess.run(["uv", "pip", "install", "--python", str(venv_py),
                         str(repo_dir)],
                        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: uv pip install of sushie failed:\n{r.stderr}")

    # cyvcf2 ships a prebuilt macOS wheel linked against an older htslib;
    # on hosts where the local htslib differs we hit
    # `symbol not found in flat namespace '_bcf_float_missing'`. Force a
    # source rebuild of cyvcf2 against the local htslib install
    # (homebrew/conda) so the dynamic symbols match.
    r = subprocess.run(
        [str(venv_py), "-c",
         "from cyvcf2 import VCF"],
        capture_output=True, text=True)
    if r.returncode != 0 and "bcf_float_missing" in r.stderr:
        print("[sushie] cyvcf2 dylib mismatch detected; rebuilding from "
              "source against local htslib ...", file=sys.stderr)
        # cyvcf2's source build calls `autoreconf -i`; ensure it's there.
        if not shutil.which("autoreconf"):
            if sys.platform == "darwin" and shutil.which("brew"):
                print("[sushie] installing autoconf/automake/libtool ...",
                      file=sys.stderr)
                subprocess.run(["brew", "install", "autoconf", "automake",
                                  "libtool"], capture_output=True, text=True)
            elif shutil.which("apt-get"):
                subprocess.run(["sudo", "apt-get", "install", "-y",
                                  "autoconf", "automake", "libtool"],
                                capture_output=True, text=True)
        env = os.environ.copy()
        # Homebrew prefix on Apple Silicon vs Intel macs.
        prefix = "/opt/homebrew" if Path("/opt/homebrew").exists() \
            else "/usr/local"
        # cyvcf2 1.x ships its own htslib by default; tell it to use the
        # system install (homebrew or conda) instead so the linked
        # symbols match the dylib actually on the host.
        env["CYVCF2_HTSLIB_MODE"] = "EXTERNAL"
        env["HTSLIB_LIBRARY_DIR"] = f"{prefix}/lib"
        env["HTSLIB_INCLUDE_DIR"] = f"{prefix}/include"
        env["CFLAGS"] = (env.get("CFLAGS", "") +
                          f" -I{prefix}/include").strip()
        env["LDFLAGS"] = (env.get("LDFLAGS", "") +
                           f" -L{prefix}/lib").strip()
        # Sushie pins pandas==1.5.0 which doesn't support numpy>=2.0
        # ('numpy.dtype size changed' ABI break). Pin numpy<2 in the
        # cyvcf2 rebuild so the venv's pandas keeps working.
        r2 = subprocess.run(
            ["uv", "pip", "install", "--python", str(venv_py),
             "--force-reinstall", "--no-binary", "cyvcf2",
             "numpy<2", "cyvcf2"],
            capture_output=True, text=True, env=env)
        if r2.returncode != 0:
            sys.exit(f"ERROR: cyvcf2 source rebuild failed:\n{r2.stderr}")
    return venv_py


def get_repo_commit(repo_dir: Path) -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo_dir,
                       capture_output=True, text=True)
    return r.stdout.strip()


def which_sushie_cli(venv_py: Optional[Path] = None) -> str:
    """Return the command (as a space-separated string for shell-equivalent
    invocation, split before subprocess.call) that launches sushie.

    Prefers (1) the `sushie` console-script in the venv's bin/ directory
    when an isolated venv was created, (2) a `sushie` on PATH, (3) the
    host `python -m sushie` as last resort (rarely works since sushie
    doesn't ship a __main__).
    """
    if venv_py is not None:
        venv_sushie = venv_py.parent / "sushie"
        if venv_sushie.exists():
            return str(venv_sushie)
    p = shutil.which("sushie")
    return p if p else f"{sys.executable} -m sushie"


# ---------------------------- Verify install ----------------------------

def verify_install(repo_dir: Path, venv_py: Optional[Path] = None) -> int:
    """Run the bundled 3-ancestry tutorial to verify the install works."""
    sentinel = CACHE_ROOT / "verify_install" / ".verified"
    data_dir = repo_dir / "data"
    out_dir = CACHE_ROOT / "verify_install"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "test_result"
    if not data_dir.exists():
        sys.exit(f"ERROR: bundled tutorial data not found at {data_dir}")
    cli = which_sushie_cli(venv_py)
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

def run_finemap(args, mode: str, repo_dir: Path,
                 venv_py: Optional[Path] = None) -> int:
    """Invoke `sushie finemap` with the user's flags. Mode is 'susie' or
    'sushie' — chosen by number of VCF args."""
    cli = which_sushie_cli(venv_py).split()
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


def _find_sushie_output(out_prefix: Path, suffix: str) -> Optional[Path]:
    """Sushie writes `<out_prefix>.sushie.<suffix>` (e.g. .sushie.cs.tsv)
    for the `region` subcommand and `<out_prefix>.<suffix>` for the
    individual-level tutorial path. Return the existing path, or None.
    """
    for candidate in (Path(f"{out_prefix}.sushie.{suffix}"),
                       Path(f"{out_prefix}.{suffix}")):
        if candidate.exists():
            return candidate
    return None


def _parse_outputs(out_prefix: Path) -> dict:
    """Parse sushie's per-CS + weights tables. Handles both naming
    conventions: the `region` subcommand produces `<prefix>.sushie.cs.tsv`
    + `<prefix>.sushie.weights.tsv`; the individual-level tutorial
    produces `<prefix>.cs.tsv` + `<prefix>.weight.tsv`."""
    summary: dict = {}
    cs_path = _find_sushie_output(out_prefix, "cs.tsv")
    if cs_path is not None:
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
            if rows and "pip_all" in header:
                pip_idx = header.index("pip_all")
                pips = []
                for r in rows:
                    try:
                        pips.append(float(r[pip_idx]))
                    except (ValueError, IndexError):
                        continue
                if pips:
                    summary["max_pip"] = round(max(pips), 4)
        except Exception as e:
            summary["cs_parse_error"] = str(e)
    weight_path = _find_sushie_output(out_prefix, "weights.tsv") or \
                   _find_sushie_output(out_prefix, "weight.tsv")
    if weight_path is not None:
        try:
            with open(weight_path) as f:
                header = f.readline().rstrip("\n").split("\t")
                # The region subcommand uses `sushie_pip_all`; the
                # tutorial uses `PIP`. Prefer whichever exists.
                for pip_name in ("sushie_pip_all", "PIP", "pip_all"):
                    if pip_name in header:
                        pip_idx = header.index(pip_name)
                        break
                else:
                    pip_idx = None
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
    corr_path = _find_sushie_output(out_prefix, "corr.tsv")
    if corr_path is not None:
        summary["cross_ancestry_correlation_file"] = str(corr_path)
    return summary


# ---------------------------- htslib tool install + tabix slicing ----------------------------

_HTSLIB_BINS = ("tabix", "bgzip", "bcftools")
_FINEMAP_REFCACHE = Path.home() / ".cache" / "finemap"

# Per-chromosome 1000G phase3 VCF URL pattern (build hg19).
_1000G_HG19_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/"
                    "20130502/ALL.chr{chrom}.phase3_shapeit2_mvncall_"
                    "integrated_v5b.20130502.genotypes.vcf.gz")
# Build GRCh38 NYGC re-call:
_1000G_HG38_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
                    "data_collections/1000_genomes_project/release/"
                    "20190312_biallelic_SNV_and_INDEL/"
                    "ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_"
                    "27022019.GRCh38.phased.vcf.gz")
# Sample → population integrated panel file.
_1000G_PANEL_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/"
                     "20130502/integrated_call_samples_v3.20130502."
                     "ALL.panel")

_POP_TO_SUPERPOP = {
    "eur": "EUR", "eas": "EAS", "afr": "AFR", "sas": "SAS", "amr": "AMR",
}


def ensure_htslib_tools() -> dict:
    """Make sure `tabix`, `bgzip`, `bcftools` are on PATH. If missing, try
    a platform-appropriate install (brew on macOS, apt on Debian/Ubuntu,
    conda if no system path works). Returns {bin: path} or sys.exits.
    """
    paths = {b: shutil.which(b) for b in _HTSLIB_BINS}
    missing = [b for b, p in paths.items() if not p]
    if not missing:
        return paths

    print(f"[finemap] installing missing htslib tools {missing} ...",
          file=sys.stderr)
    # Install htslib/bcftools + autoconf/automake/libtool — the latter
    # three are needed when cyvcf2 has to rebuild from source against
    # the local htslib (its CMake calls `autoreconf -i`).
    installers = []
    if sys.platform == "darwin" and shutil.which("brew"):
        installers.append(["brew", "install", "htslib", "bcftools",
                            "autoconf", "automake", "libtool"])
    if shutil.which("apt-get"):
        installers.append(["sudo", "apt-get", "install", "-y",
                            "tabix", "bcftools", "autoconf", "automake",
                            "libtool"])
    if shutil.which("conda"):
        installers.append(["conda", "install", "-y", "-c", "bioconda",
                            "htslib", "bcftools"])
    if shutil.which("mamba"):
        installers.append(["mamba", "install", "-y", "-c", "bioconda",
                            "htslib", "bcftools"])

    for cmd in installers:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                print(f"[finemap] installed via: {' '.join(cmd[:3])}",
                      file=sys.stderr)
                break
        except FileNotFoundError:
            continue

    paths = {b: shutil.which(b) for b in _HTSLIB_BINS}
    still_missing = [b for b, p in paths.items() if not p]
    if still_missing:
        sys.exit(f"ERROR: could not install {still_missing}. Install "
                 f"manually: `brew install htslib bcftools` (macOS) / "
                 f"`apt-get install tabix bcftools` (Debian) / "
                 f"`conda install -c bioconda htslib bcftools`.")
    return paths


def _fetch_1000g_panel(cache_root: Path) -> Path:
    """Cache the 1000G integrated panel file (sample → population)."""
    p = cache_root / "1000g_integrated_panel.tsv"
    if p.exists() and p.stat().st_size > 0:
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    print(f"[finemap] downloading 1000G integrated panel -> {p}",
          file=sys.stderr)
    urllib.request.urlretrieve(_1000G_PANEL_URL, p)
    return p


def _samples_for_superpop(panel_file: Path, superpop: str) -> list[str]:
    """Return the list of 1000G sample IDs in a given super-population."""
    samples = []
    with open(panel_file) as f:
        header = f.readline().rstrip("\n").split()
        try:
            i_sample = header.index("sample")
            i_super = header.index("super_pop")
        except ValueError:
            return []
        for line in f:
            r = line.rstrip("\n").split()
            if len(r) > i_super and r[i_super] == superpop:
                samples.append(r[i_sample])
    return samples


def ensure_1000g_slice(chrom: int, start: int, end: int,
                        population: str, build: str = "hg19",
                        cache_root: Path = _FINEMAP_REFCACHE) -> Path:
    """Tabix-slice the 1000G phase3 VCF to (chrom:start-end), subset to a
    super-population, bgzip+index, and cache. Returns the local cached path.

    Args:
        chrom: 1-22 (numeric)
        start, end: bp coordinates on `build`
        population: eur, eas, afr, sas, amr (case-insensitive)
        build: hg19 (default) or hg38

    Reuses the cached file when present. Total fetch on first run is
    typically a few MB for a 1 Mb window.
    """
    pop = population.lower()
    superpop = _POP_TO_SUPERPOP.get(pop)
    if not superpop:
        sys.exit(f"REFUSED: unknown population {population!r}; "
                 f"allowed: {list(_POP_TO_SUPERPOP)}")
    if build not in ("hg19", "hg38"):
        sys.exit(f"REFUSED: unknown build {build!r}; allowed: hg19, hg38")

    ensure_htslib_tools()
    cache_root = cache_root / f"1000g_{build}_{pop}"
    cache_root.mkdir(parents=True, exist_ok=True)
    out_vcf = cache_root / f"chr{chrom}_{start}_{end}.vcf.gz"
    if out_vcf.exists() and out_vcf.stat().st_size > 0:
        idx = Path(str(out_vcf) + ".tbi")
        if idx.exists():
            return out_vcf

    url_template = _1000G_HG19_URL if build == "hg19" else _1000G_HG38_URL
    url = url_template.format(chrom=chrom)
    region = f"{chrom}:{start}-{end}"

    # 1) Cache the sample → population panel.
    panel = _fetch_1000g_panel(_FINEMAP_REFCACHE)
    samples = _samples_for_superpop(panel, superpop)
    if not samples:
        sys.exit(f"REFUSED: no {superpop} samples found in 1000G panel "
                 f"at {panel}")
    keep_file = cache_root / f"{pop}.samples.txt"
    keep_file.write_text("\n".join(samples) + "\n")

    # 2) tabix-slice the remote VCF to the window, pipe into bcftools to
    # subset samples, write bgzipped output, then tabix-index.
    print(f"[finemap] tabix-slicing {url} @ {region} ({len(samples)} "
          f"{superpop} samples) -> {out_vcf}", file=sys.stderr)
    tmp_raw = cache_root / f"chr{chrom}_{start}_{end}.raw.vcf.gz"
    # Stream tabix output through bgzip directly so we never materialise the
    # full-population intermediate as plain VCF.
    with open(tmp_raw, "wb") as out:
        p_tab = subprocess.Popen(["tabix", "-h", url, region],
                                  stdout=subprocess.PIPE)
        p_bz = subprocess.Popen(["bgzip", "-c"],
                                 stdin=p_tab.stdout, stdout=out)
        p_tab.stdout.close()
        p_bz.communicate()
        p_tab.wait()
        if p_tab.returncode != 0 or p_bz.returncode != 0:
            sys.exit(f"ERROR: tabix/bgzip pipe failed for {region} from "
                     f"{url}")
    # bcftools view -S to subset samples + annotate ID = chr:pos. The
    # GWAS converter writes IDs in the same `chr:pos` form, so the join
    # is position-only — sushie's internal logic does allele alignment
    # from the genotype dosages vs the GWAS effect/other alleles.
    # Filter to biallelic SNPs (-m2 -M2 -v snps) to drop multi-allelic
    # sites that would create duplicate IDs.
    tmp_sub = out_vcf.with_suffix(".sub.vcf.gz")
    rc = subprocess.call(["bcftools", "view",
                           "-S", str(keep_file),
                           "-m2", "-M2", "-v", "snps",
                           "-Oz", "-o", str(tmp_sub),
                           str(tmp_raw)])
    if rc != 0:
        sys.exit(f"ERROR: bcftools view failed (rc={rc})")
    tmp_raw.unlink(missing_ok=True)
    rc = subprocess.call([
        "bcftools", "annotate",
        "--set-id", "%CHROM:%POS",
        "-Oz", "-o", str(out_vcf), str(tmp_sub)])
    if rc != 0:
        sys.exit(f"ERROR: bcftools annotate failed (rc={rc})")
    tmp_sub.unlink(missing_ok=True)

    rc = subprocess.call(["tabix", "-p", "vcf", str(out_vcf)])
    if rc != 0:
        sys.exit(f"ERROR: tabix index failed for {out_vcf}")
    return out_vcf


# ---------------------------- GWAS region mode (sumstats + ref VCF/PLINK) ----------------------------

# sushie's `finemap --summary` mode reads a sumstats TSV with this exact
# column order; we auto-convert from GWAS-Catalog harmonised or LDSC
# munged input.
SUSHIE_GWAS_HEADER = ["chrom", "snp", "pos", "a1", "a0", "z"]

_NA_TOKENS = {"", "NA", "nan", "NaN", "N/A", ".", "null", "None"}


def _open_text(path: Path):
    import gzip
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def detect_and_convert_to_sushie_gwas(input_path: Path,
                                       chrom: int, start: int, end: int,
                                       out_path: Path,
                                       extras_path: Optional[Path] = None
                                       ) -> int:
    """Read a GWAS sumstats file (GWAS Catalog harmonised TSV, LDSC munged
    `.sumstats.gz`, or any TSV with chrom/snp/pos/effect_allele/other_allele
    + signed effect or z columns), filter to the (chrom, start, end) window,
    and write a sushie-compatible TSV with header
    `chrom snp pos a1 a0 z`.

    Returns the number of rows written.
    """
    with _open_text(input_path) as f:
        header = f.readline().rstrip("\n").split("\t")
    cols = {c: i for i, c in enumerate(header)}
    cols_l = {c.lower(): i for i, c in enumerate(header)}

    def col(*candidates):
        for c in candidates:
            if c in cols:
                return cols[c]
            if c.lower() in cols_l:
                return cols_l[c.lower()]
        return None

    # Prefer harmonised columns; alleles get uppercased on the way out.
    chr_i = col("CHR", "chr", "chromosome", "hm_chrom")
    pos_i = col("POS", "pos", "base_pair_location", "hm_pos", "bp")
    snp_i = col("SNP", "snp", "hm_rsid", "rsid", "variant_id", "snpid")
    a1_i  = col("A1", "a1", "hm_effect_allele", "effect_allele")
    a0_i  = col("A2", "a0", "a2", "hm_other_allele", "other_allele")
    beta_i = col("BETA", "beta", "hm_beta", "b")
    se_i   = col("SE", "se", "standard_error", "stderr")
    z_i    = col("Z", "z", "zscore")
    or_i   = col("OR", "hm_odds_ratio", "odds_ratio")
    p_i    = col("P", "p", "p_value", "pval", "Pvalue", "PValue")

    missing = [name for name, idx in
                [("chrom", chr_i), ("pos", pos_i), ("snp", snp_i),
                 ("a1", a1_i), ("a0", a0_i)] if idx is None]
    if missing:
        sys.exit(f"REFUSED: input header missing required columns "
                 f"{missing}; header was {header[:25]}")
    has_signed = z_i is not None or (beta_i is not None and se_i is not None) \
                  or (or_i is not None and se_i is not None)
    if not has_signed:
        sys.exit("REFUSED: cannot derive a signed effect from the input "
                 "(need Z, or BETA+SE, or OR+SE)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    # Sidecar `<out>.gwas_extras.tsv` carries the per-SNP rsid + beta + se
    # + p mapped by chr:pos, so the post-run CS table can join back to
    # rsids and effect-size details (sushie's CS output uses chr:pos as
    # the SNP id and doesn't preserve rsids).
    extras_handle = None
    if extras_path is not None:
        extras_path.parent.mkdir(parents=True, exist_ok=True)
        extras_handle = open(extras_path, "w")
        extras_handle.write("chrpos\trsid\tbeta\tse\tp\n")
    with _open_text(input_path) as f, open(out_path, "w") as g:
        f.readline()  # skip header
        g.write("\t".join(SUSHIE_GWAS_HEADER) + "\n")
        for line in f:
            r = line.rstrip("\n").split("\t")
            if len(r) < max(i for i in (chr_i, pos_i, snp_i, a1_i, a0_i,
                                          beta_i, se_i, z_i, or_i)
                              if i is not None) + 1:
                continue
            # Normalise chrom to int; skip X/Y/MT (sushie only does 1-22).
            chrom_s = str(r[chr_i]).lstrip("chr").strip()
            if not chrom_s.isdigit():
                continue
            try:
                pos_int = int(r[pos_i])
            except ValueError:
                continue
            chrom_i = int(chrom_s)
            if chrom_i != int(chrom) or pos_int < int(start) or pos_int > int(end):
                continue
            snp = r[snp_i]
            a1 = r[a1_i].upper(); a0 = r[a0_i].upper()
            if snp in _NA_TOKENS or a1 in _NA_TOKENS or a0 in _NA_TOKENS:
                continue
            # signed effect -> z
            try:
                if z_i is not None and r[z_i] not in _NA_TOKENS:
                    z = float(r[z_i])
                elif beta_i is not None and se_i is not None:
                    b = float(r[beta_i]); se = float(r[se_i])
                    if se == 0:
                        continue
                    z = b / se
                elif or_i is not None and se_i is not None:
                    import math
                    b = math.log(float(r[or_i])); se = float(r[se_i])
                    if se == 0:
                        continue
                    z = b / se
                else:
                    continue
            except (ValueError, IndexError):
                continue
            # Standardise SNP ID to plain `chr:pos`. The 1000G phase3 VCFs
            # ship with ID=. so rsid joins fail; using alleles in the ID
            # (chr:pos:A1:A0) is too strict — strand flips, allele
            # ordering, and multi-allelic encodings drop most SNPs at the
            # join. Position-only IDs maximise the SNP overlap; sushie's
            # internal logic then handles allele alignment.
            snp_id = f"{chrom_i}:{pos_int}"
            g.write(f"{chrom_i}\t{snp_id}\t{pos_int}\t{a1}\t{a0}\t{z:.6g}\n")
            if extras_handle is not None:
                # Capture beta + se + p if present, else NA — used by the
                # post-run summary table to reconstruct effect-size detail.
                def _safe(idx):
                    if idx is None:
                        return "NA"
                    try:
                        v = r[idx]
                        return v if v not in _NA_TOKENS else "NA"
                    except IndexError:
                        return "NA"
                p_str = _safe(p_i)
                if beta_i is not None and se_i is not None:
                    b_str = _safe(beta_i); se_str = _safe(se_i)
                elif or_i is not None and se_i is not None:
                    try:
                        b_str = f"{math.log(float(r[or_i])):.6g}"
                    except (ValueError, IndexError):
                        b_str = "NA"
                    se_str = _safe(se_i)
                else:
                    b_str = "NA"; se_str = _safe(se_i)
                extras_handle.write(f"{snp_id}\t{snp}\t{b_str}\t{se_str}\t{p_str}\n")
            n_written += 1
    if extras_handle is not None:
        extras_handle.close()
    return n_written


def run_region(args, repo_dir: Path,
                venv_py: Optional[Path] = None) -> int:
    """Single-locus fine-mapping with sushie's `finemap --summary`.
    Sushie computes LD internally from --vcf / --plink / --bgen reference
    genotypes (no precomputed LD matrix required)."""
    sushie_gwas = Path(str(args.out) + ".gwas.tsv")
    extras_sidecar = Path(str(args.out) + ".gwas_extras.tsv")
    n = detect_and_convert_to_sushie_gwas(
        args.gwas_sumstats, args.chrom, args.start, args.end, sushie_gwas,
        extras_path=extras_sidecar)
    print(f"[finemap region] converted GWAS -> {sushie_gwas} ({n} rows "
          f"in chr{args.chrom}:{args.start}-{args.end})", file=sys.stderr)
    if n < 50:
        sys.exit(f"REFUSED: only {n} GWAS rows in the requested window; "
                 f"need at least 50 SNPs to fine-map")

    cli = which_sushie_cli(venv_py).split()
    cmd = cli + ["finemap",
                  "--summary",
                  "--gwas", str(sushie_gwas),
                  "--gwas-header", *SUSHIE_GWAS_HEADER,
                  "--sample-size", str(args.N),
                  "--chrom", str(args.chrom),
                  "--start", str(args.start),
                  "--end", str(args.end),
                  "--output", str(args.out)]
    # LD source resolution (in priority order):
    #   1. Explicit --ref-vcf / --ref-plink / --ref-bgen / --ld
    #   2. --ref-1000g <pop>: tabix-slice from 1000G FTP + subset samples
    ref_vcf_paths = list(args.ref_vcf or [])
    if args.ref_1000g and not ref_vcf_paths and not args.ref_plink and \
            not args.ref_bgen and not args.ld:
        sliced = ensure_1000g_slice(args.chrom, args.start, args.end,
                                      population=args.ref_1000g,
                                      build=args.ref_1000g_build)
        ref_vcf_paths = [sliced]
        print(f"[finemap region] 1000G {args.ref_1000g.upper()} "
              f"slice -> {sliced}", file=sys.stderr)
    if ref_vcf_paths:
        cmd += ["--vcf", *[str(p) for p in ref_vcf_paths]]
    elif args.ref_plink:
        cmd += ["--plink", *[str(p) for p in args.ref_plink]]
    elif args.ref_bgen:
        cmd += ["--bgen", *[str(p) for p in args.ref_bgen]]
    elif args.ld:
        cmd += ["--ld", *[str(p) for p in args.ld]]
    else:
        sys.exit("REFUSED: provide one of --ref-vcf / --ref-plink / "
                 "--ref-bgen / --ld / --ref-1000g for the LD source")
    if args.extra:
        # Drop the literal `--` separator that argparse REMAINDER captures
        # but the underlying sushie CLI doesn't accept.
        cmd += [a for a in args.extra if a != "--"]

    print(f"[finemap region] {' '.join(cmd)}", file=sys.stderr)
    rc = subprocess.call(cmd)

    # OKG provenance (single-ancestry SuSiE-RSS path).
    nids = list(OKG_NODES_SUSIE.values())
    found = resolve_okg(args.okg_repo, nids)
    okg_node_ids = {k: v for k, v in OKG_NODES_SUSIE.items() if v in found}
    if args.okg_ld_panel_id:
        okg_node_ids["ld_panel"] = args.okg_ld_panel_id
    if args.okg_dataset_id:
        okg_node_ids["dataset"] = args.okg_dataset_id

    summary = _parse_outputs(args.out)
    summary["n_snps_in_window"] = n
    manifest = {
        "subcommand": "region",
        "output_prefix": str(args.out),
        "sushie_repo": args.repo_url,
        "sushie_commit": get_repo_commit(repo_dir),
        "n_ancestries": 1,
        "locus": {"chrom": int(args.chrom),
                   "start": int(args.start),
                   "end": int(args.end)},
        "inputs": {
            "gwas_sumstats": str(args.gwas_sumstats),
            "ref_vcf":   [str(p) for p in (args.ref_vcf or [])],
            "ref_plink": [str(p) for p in (args.ref_plink or [])],
            "ref_bgen":  [str(p) for p in (args.ref_bgen or [])],
            "ld":        [str(p) for p in (args.ld or [])],
        },
        "sample_size": args.N,
        "summary": summary,
        "okg_node_ids": okg_node_ids,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{args.out}.finemap.json")
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    if rc == 0 and summary:
        print(f"summary: {summary.get('n_cs', '?')} CSs, "
              f"mean size {summary.get('mean_cs_size', '?')}, "
              f"max PIP {summary.get('max_pip', '?')}")
        # Print the rsID-enriched credible-set table, optionally with
        # gnomAD + GTEx annotations appended.
        annot_sources: set = set()
        if args.annotate in ("gnomad", "both"):
            annot_sources.add("gnomad")
        if args.annotate in ("gtex", "both"):
            annot_sources.add("gtex")
        _print_cs_summary_table(Path(args.out), extras_sidecar,
                                  annot_sources=annot_sources)
    return rc


_ANNOT_CACHE = Path.home() / ".cache" / "variant-annotate"


def _load_variant_annotate():
    """Import variant-annotate/scripts/annotate.py as a module. The
    variant-annotate skill lives as a sibling directory under
    statgen-skills/. Returns the module, or None if not findable."""
    # finemap.py is at: statgen-skills/finemap/scripts/finemap.py
    # variant-annotate annotate.py is at:
    #         statgen-skills/variant-annotate/scripts/annotate.py
    here = Path(__file__).resolve().parent  # statgen-skills/finemap/scripts
    sibling = here.parent.parent / "variant-annotate" / "scripts"
    if not (sibling / "annotate.py").exists():
        return None
    if str(sibling) not in sys.path:
        sys.path.insert(0, str(sibling))
    try:
        import annotate as va_mod  # noqa
        return va_mod
    except Exception as e:
        sys.stderr.write(f"[finemap annotate] failed to import "
                          f"variant-annotate: {type(e).__name__}: {e}\n")
        return None


def _annotate_variants(rsids: list, sources: set) -> dict:
    """Fetch gnomAD and/or GTEx annotations for rsIDs by delegating to
    the sibling `variant-annotate` skill. Returns
    {rsid: {gnomad: dict, gtex: dict}}, empty dict if the skill isn't
    available. Per-rsID JSON is cached by variant-annotate under
    ~/.cache/variant-annotate/."""
    va = _load_variant_annotate()
    if va is None:
        sys.stderr.write("[finemap annotate] variant-annotate skill not "
                          "found alongside finemap; skipping annotations\n")
        return {}
    records = va.annotate_rsids(list(rsids), set(sources), _ANNOT_CACHE)
    out: dict = {}
    for rec in records:
        rsid = rec.get("rsid")
        bundle: dict = {}
        if "gnomad" in sources and rec.get("gnomad"):
            bundle["gnomad"] = rec["gnomad"]
        if "gtex" in sources and rec.get("gtex"):
            bundle["gtex"] = rec["gtex"]
        if rsid and bundle:
            out[rsid] = bundle
    return out


def _print_cs_summary_table(out_prefix: Path,
                             extras_sidecar: Path,
                             annot_sources: Optional[set] = None) -> None:
    """Print a human-readable rsID-enriched credible-set table to stdout.

    Reads sushie's `<prefix>.sushie.cs.tsv` (which carries CSIndex, snp,
    pos, a0, a1, pip_all) and joins by chr:pos with the `<prefix>.
    gwas_extras.tsv` sidecar (chrpos, rsid, beta, se, p) the converter
    wrote. Output: one row per CS member, sorted by CS index then by
    descending pip_all within each CS."""
    cs_path = _find_sushie_output(out_prefix, "cs.tsv")
    if cs_path is None or not extras_sidecar.exists():
        return
    extras = {}
    with open(extras_sidecar) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_chrpos = header.index("chrpos")
        idx_rsid = header.index("rsid")
        idx_beta = header.index("beta")
        idx_se = header.index("se")
        idx_p = header.index("p")
        for line in f:
            r = line.rstrip("\n").split("\t")
            extras[r[idx_chrpos]] = (r[idx_rsid], r[idx_beta], r[idx_se],
                                       r[idx_p])
    rows = []
    with open(cs_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_cs = header.index("CSIndex")
            i_snp = header.index("snp")
            i_pos = header.index("pos")
            i_a0 = header.index("a0")
            i_a1 = header.index("a1")
            i_pip = header.index("pip_all")
        except ValueError:
            return
        for line in f:
            r = line.rstrip("\n").split("\t")
            cs = int(r[i_cs])
            try:
                pip = float(r[i_pip])
            except ValueError:
                pip = 0.0
            chrpos = r[i_snp]
            ex = extras.get(chrpos, ("NA", "NA", "NA", "NA"))
            rows.append((cs, -pip, chrpos, r[i_pos], r[i_a0], r[i_a1],
                          pip, *ex))
    if not rows:
        return
    rows.sort(key=lambda x: (x[0], x[1]))

    # Optional: fetch gnomAD + GTEx annotations for every CS-member rsid.
    annots = {}
    if annot_sources:
        rsids = sorted({r[7] for r in rows
                        if r[7] and r[7] != "NA" and r[7].startswith("rs")})
        if rsids:
            print(f"[finemap annotate] querying "
                  f"{','.join(sorted(annot_sources))} for {len(rsids)} "
                  f"variants (cached to {_ANNOT_CACHE}) ...",
                  file=sys.stderr)
            annots = _annotate_variants(rsids, annot_sources)

    # Render markdown-style table.
    print()
    print("Credible-set members (sorted by CS, descending PIP within CS):")
    print()
    base_cols = ["CS", "PIP", "rsID", "chr:pos", "Eff", "NonEff",
                  "Beta", "SE", "P"]
    extra_cols = []
    if "gnomad" in (annot_sources or set()):
        extra_cols += ["Csq", "Gene", "AF_nfe"]
    if "gtex" in (annot_sources or set()):
        extra_cols += ["top eQTL", "top sQTL"]
    cols = base_cols + extra_cols
    print("| " + " | ".join(cols) + " |")
    print("|" + "|".join("---" for _ in cols) + "|")
    for cs, _, chrpos, pos, a0, a1, pip, rsid, beta, se, p in rows:
        row = [str(cs), f"{pip:.4f}", rsid, chrpos, a1, a0,
                beta, se, p]
        rec = annots.get(rsid, {})
        if "gnomad" in (annot_sources or set()):
            g = rec.get("gnomad") or {}
            csq = (g.get("consequence") or "").replace("_variant", "")
            gene = g.get("gene") or ""
            af_nfe = g.get("af", {}).get("nfe")
            af_str = f"{af_nfe:.4g}" if isinstance(af_nfe, (int, float)) \
                else ""
            row += [csq, gene, af_str]
        if "gtex" in (annot_sources or set()):
            t = rec.get("gtex") or {}
            eq = t.get("top_eqtl") or {}
            sq = t.get("top_sqtl") or {}
            eq_str = (f"{eq.get('gene')}@{eq.get('tissue')} "
                       f"p={eq.get('p'):.1e}") if eq.get("gene") else ""
            sq_str = (f"{sq.get('gene')}@{sq.get('tissue')} "
                       f"p={sq.get('p'):.1e}") if sq.get("gene") else ""
            row += [eq_str, sq_str]
        print("| " + " | ".join(row) + " |")
    print()


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

    preg = sub.add_parser("region",
                           help="Single-locus fine-mapping from GWAS sumstats "
                                "+ a reference VCF/PLINK (sushie computes "
                                "LD internally; no precomputed matrix needed)")
    preg.add_argument("--gwas-sumstats", dest="gwas_sumstats", type=Path,
                      required=True,
                      help="GWAS sumstats file (GWAS Catalog harmonised "
                           "TSV, LDSC munged .sumstats.gz, or any TSV "
                           "with chrom/snp/pos/effect_allele/other_allele "
                           "+ z or beta/se)")
    preg.add_argument("--chrom", type=int, required=True, choices=range(1, 23),
                      help="Chromosome (1-22)")
    preg.add_argument("--start", type=int, required=True,
                      help="Window start, bp")
    preg.add_argument("--end", type=int, required=True,
                      help="Window end, bp")
    preg.add_argument("--N", type=int, required=True,
                      help="GWAS total sample size for this trait")
    preg.add_argument("--ref-vcf", dest="ref_vcf", type=Path, nargs="+",
                      help="Reference genotype VCF (sushie computes LD "
                           "from this; usually 1000G or UKB on a build "
                           "matching the GWAS)")
    preg.add_argument("--ref-plink", dest="ref_plink", type=Path, nargs="+",
                      help="Reference genotype PLINK1.9 prefix(es) "
                           "(alternative to --ref-vcf)")
    preg.add_argument("--ref-bgen", dest="ref_bgen", type=Path, nargs="+",
                      help="Reference genotype BGEN 1.3 file(s) "
                           "(alternative to --ref-vcf)")
    preg.add_argument("--ref-1000g", dest="ref_1000g", type=str,
                      choices=list(_POP_TO_SUPERPOP),
                      help="Convenience: auto-fetch the 1000G phase3 "
                           "VCF slice for the given super-population "
                           "(eur/eas/afr/sas/amr) at --chrom:--start-end, "
                           "subset to that population, and use it as the "
                           "LD reference. Installs htslib + bcftools on "
                           "first use; caches the slice under "
                           "~/.cache/finemap/1000g_<build>_<pop>/.")
    preg.add_argument("--ref-1000g-build", dest="ref_1000g_build",
                      type=str, choices=["hg19", "hg38"], default="hg19",
                      help="Build for --ref-1000g (default: hg19)")
    preg.add_argument("--ld", type=Path, nargs="+",
                      help="Pre-computed LD matrix (tsv/tsv.gz) — bypasses "
                           "the reference-genotype LD computation")
    preg.add_argument("--out", type=Path, required=True,
                      help="Output prefix")
    preg.add_argument("--okg-dataset-id", dest="okg_dataset_id", type=str,
                      help="OKG dataset_metadata node ID for the GWAS "
                           "(recorded in the manifest's okg_node_ids)")
    preg.add_argument("--okg-ld-panel-id", dest="okg_ld_panel_id", type=str,
                      help="OKG ld_panel node ID for the reference "
                           "genotypes (recorded in the manifest)")
    preg.add_argument("--annotate", type=str, default="both",
                      choices=["none", "gnomad", "gtex", "both"],
                      help="Post-fine-mapping variant annotation, delegated "
                           "to the sibling `variant-annotate` skill. "
                           "`gnomad` adds most-severe consequence + gene + "
                           "NFE allele frequency. `gtex` adds top eQTL + "
                           "sQTL across tissues. `both` (default) adds both "
                           "column groups to the CS summary table. Network "
                           "calls; responses cached at "
                           "~/.cache/variant-annotate/.")
    preg.add_argument("--repo-url", type=str, default=DEFAULT_REPO)
    preg.add_argument("--sushie-commit", dest="commit", type=str, default=None)
    preg.add_argument("--repo-cache", type=Path, default=CACHE_ROOT)
    preg.add_argument("--refresh", action="store_true")
    preg.add_argument("--re-verify", action="store_true")
    preg.add_argument("--okg-repo", type=Path,
                      default=Path(os.environ["OKG_REPO"])
                              if os.environ.get("OKG_REPO") else None)
    preg.add_argument("extra", nargs=argparse.REMAINDER,
                      help="Extra flags forwarded to `sushie finemap`")

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
        venv_py = ensure_sushie_installed(repo_dir)
        return verify_install(repo_dir, venv_py=venv_py)

    if args.subcmd not in ("susie", "sushie", "sumstats", "region"):
        p.print_help(sys.stderr)
        return 2

    # Re-parse for the subcommand's flags.
    args = p.parse_args()

    # Common setup.
    repo_dir = ensure_repo(args.repo_url, args.commit,
                            Path(args.repo_cache), refresh=args.refresh)
    venv_py = ensure_sushie_installed(repo_dir)

    # First-run verify-install gate.
    if args.refresh or args.re_verify or not is_verified():
        rc = verify_install(repo_dir, venv_py=venv_py)
        if rc != 0:
            sys.exit("verify-install failed; refusing to proceed with user run")

    if args.subcmd == "sumstats":
        return run_sumstats(args, repo_dir)
    if args.subcmd == "region":
        return run_region(args, repo_dir, venv_py=venv_py)

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

    return run_finemap(args, args.subcmd, repo_dir, venv_py=venv_py)


if __name__ == "__main__":
    sys.exit(main())
