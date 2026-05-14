#!/usr/bin/env python3
"""Polygenic risk score (PRS) skill orchestrator.

Currently implements one PRS method behind a `--method` flag:

    prs.py --method sbayesrc --gwas-sumstats <path> --ancestry eur \\
           --out <prefix> [--okg-dataset-id dataset:...]

OKG resolution:
  - --ancestry {eur,eas,afr}     -> ld_panel:sbayesrc_hm3_<anc>
  - --okg-ld-panel-id <node_id>  -> explicit override
  - --okg-dataset-id <node_id>   -> reads n_cases + n_controls (or n_samples)
                                    to populate `N` in the COJO file.

The skill writes a `.prs.json` sidecar that cites the method/software/paper/
ld_panel/dataset OKG node IDs plus the SHA-256 of the LD reference.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


CACHE_ROOT = Path.home() / ".cache" / "sbayesrc"

LD_PANEL_BY_ANC = {
    "eur": "ld_panel:sbayesrc_hm3_eur",
    "eas": "ld_panel:sbayesrc_hm3_eas",
    "afr": "ld_panel:sbayesrc_hm3_afr",
}

OKG_NODES_SBAYESRC = {
    "method": "method:sbayesrc",
    "software": "software:sbayesrc_r",
    "paper": "paper:zheng_2024_sbayesrc",
}

# Columns SBayesRC expects in COJO format.
COJO_COLS = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]

# Tokens that mean "NA" in the harmonised input.
_NA_TOKENS = {"", "NA", "nan", "NaN", "N/A", ".", "null", "None"}


# --------------------------- OKG resolution ---------------------------

def _mcp_call(okg_repo: Path, method: str, arguments: dict) -> Optional[dict]:
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
                         "clientInfo": {"name": "prs-skill", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": method, "arguments": arguments}})
        resp = read()
    except Exception:
        return None
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    return resp.get("result", {}).get("structuredContent") or None


def resolve_ld_panel(okg_repo: Optional[Path],
                     ancestry: Optional[str],
                     panel_id: Optional[str]) -> dict:
    """Return {node_id, source_url, local_path_hint, genome_build, ancestry}
    for the LD panel. Raises SystemExit on REFUSED conditions."""
    if not panel_id:
        if not ancestry:
            sys.exit("REFUSED: provide either --ancestry {eur,eas,afr} or "
                     "--okg-ld-panel-id <node_id>")
        panel_id = LD_PANEL_BY_ANC.get(ancestry.lower())
        if not panel_id:
            sys.exit(f"REFUSED: unknown ancestry {ancestry!r}. "
                     f"Allowed: {list(LD_PANEL_BY_ANC)}")
    if okg_repo is None:
        sys.exit("REFUSED: $OKG_REPO not set; can't resolve LD panel metadata "
                 "from the graph. Pass --okg-repo or set the env var.")
    sc = _mcp_call(okg_repo, "get_node", {"node_id": panel_id})
    if not sc or sc.get("node_id") != panel_id:
        sys.exit(f"REFUSED: OKG has no node {panel_id} at the current "
                 f"generation. Add it via an OpenSpec change first.")
    a = sc.get("attrs") or {}
    return {
        "node_id": panel_id,
        "source_url": a.get("source_url"),
        "local_path_hint": a.get("local_path_hint"),
        "genome_build": a.get("genome_build"),
        "ancestry": a.get("ancestry_scope"),
        "name": a.get("name"),
    }


def resolve_dataset_n(okg_repo: Optional[Path],
                      dataset_id: Optional[str]) -> Optional[int]:
    if okg_repo is None or not dataset_id:
        return None
    sc = _mcp_call(okg_repo, "get_node", {"node_id": dataset_id})
    if not sc or sc.get("node_id") != dataset_id:
        return None
    a = sc.get("attrs") or {}
    try:
        if a.get("n_cases") is not None and a.get("n_controls") is not None:
            return int(a["n_cases"]) + int(a["n_controls"])
        if a.get("n_samples") is not None:
            return int(a["n_samples"])
    except (TypeError, ValueError):
        return None
    return None


# --------------------------- LD download / extract ---------------------------

def _looks_like_sbayesrc_ld_dir(d: Path) -> bool:
    """SBayesRC's LD-eigendecomp folder always contains a `snp.info` file
    plus per-block `block*.eigen.bin` artifacts. Use this as the recogniser
    so we don't return an outer wrapper dir by mistake."""
    if not d.is_dir():
        return False
    return (d / "snp.info").exists() or any(d.glob("*.eigen.bin"))


def _descend_to_ld_dir(d: Path) -> Path:
    """If `d` doesn't carry the SBayesRC LD layout itself but has exactly
    one nested directory that does, return the nested dir. Otherwise
    return `d` unchanged."""
    if _looks_like_sbayesrc_ld_dir(d):
        return d
    children = [p for p in d.iterdir() if p.is_dir()]
    if len(children) == 1 and _looks_like_sbayesrc_ld_dir(children[0]):
        return children[0]
    for c in children:
        if _looks_like_sbayesrc_ld_dir(c):
            return c
    return d


def ensure_ld_dir(panel: dict, cache_root: Path) -> Path:
    """Download (if needed) + unzip the LD eigendecomposition. Returns the
    path to the unzipped folder that actually carries SBayesRC's
    `snp.info` + `block*.eigen.bin` files (descending past wrapper dirs)."""
    hint = panel.get("local_path_hint")
    if hint:
        d = Path(os.path.expanduser(hint))
        if d.exists() and any(d.iterdir()):
            return _descend_to_ld_dir(d)
    url = panel.get("source_url")
    if not url:
        sys.exit(f"REFUSED: OKG node {panel['node_id']} has no source_url; "
                 f"cannot download LD reference.")
    cache_root.mkdir(parents=True, exist_ok=True)
    zname = Path(url).name
    zpath = cache_root / zname
    if not zpath.exists() or zpath.stat().st_size == 0:
        print(f"[prs] downloading {url} -> {zpath} "
              f"(this may take a while; ~3-5 GB)", file=sys.stderr)
        urllib.request.urlretrieve(url, zpath)
    extract_root = cache_root / zname.replace(".zip", "")
    if not extract_root.exists() or not any(extract_root.iterdir()):
        print(f"[prs] extracting {zpath} -> {extract_root}", file=sys.stderr)
        extract_root.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zpath, "r") as z:
            z.extractall(extract_root)
    # The zip often unpacks into a same-named nested folder; descend to the
    # dir that actually carries SBayesRC's snp.info + block*.eigen.bin.
    return _descend_to_ld_dir(extract_root)


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------- COJO conversion ---------------------------

def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _load_ld_snp_info(ld_dir: Path) -> dict:
    """Return {rsid: (A1, A2, A1Freq)} from the LD reference's snp.info.
    Used to fill NA freq + sanity-check alleles when the input GWAS
    doesn't carry a usable effect_allele_frequency column."""
    info_path = ld_dir / "snp.info"
    if not info_path.exists():
        return {}
    out = {}
    with open(info_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_id = header.index("ID")
            i_a1 = header.index("A1")
            i_a2 = header.index("A2")
            i_frq = header.index("A1Freq")
        except ValueError:
            return {}
        for line in f:
            r = line.rstrip("\n").split("\t")
            if len(r) > i_frq:
                try:
                    out[r[i_id]] = (r[i_a1].upper(), r[i_a2].upper(),
                                      float(r[i_frq]))
                except ValueError:
                    continue
    return out


def detect_format_and_convert_to_cojo(input_path: Path,
                                        N_const: Optional[int],
                                        out_cojo: Path,
                                        ld_dir: Optional[Path] = None) -> int:
    """Detect input format (LDSC munged, GWAS-Catalog harmonised hm_*, or
    bare GWAS-SSF) and emit a COJO TSV with header `SNP A1 A2 freq b se p N`.
    Returns the number of rows written."""
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

    # Prefer the GWAS-Catalog `hm_*` (harmonised) columns over their
    # non-hm siblings because the bare columns carry lowercase alleles
    # in some deposits, which break downstream allele-match checks.
    snp_i = col("SNP", "snp", "hm_rsid", "rsid", "variant_id", "snpid")
    a1_i  = col("A1", "a1", "hm_effect_allele", "effect_allele")
    a2_i  = col("A2", "a2", "hm_other_allele", "other_allele")
    p_i   = col("P", "p", "p_value", "pval")
    beta_i = col("BETA", "beta", "hm_beta", "b")
    se_i   = col("SE", "se", "standard_error", "stderr")
    z_i    = col("Z", "z", "zscore")
    or_i   = col("OR", "hm_odds_ratio", "odds_ratio")
    frq_i  = col("FRQ", "freq", "hm_effect_allele_frequency",
                  "effect_allele_frequency", "EAF")
    n_i    = col("N", "n", "n_total", "sample_size")
    n_cas_i = col("N_CAS", "n_cases", "ncas")
    n_con_i = col("N_CON", "n_controls", "ncon")

    if snp_i is None or a1_i is None or a2_i is None or p_i is None:
        sys.exit(f"REFUSED: cannot map COJO columns from header {header!r}. "
                 f"Need at least SNP/rsid + A1/effect_allele + "
                 f"A2/other_allele + p_value.")
    # Need a signed effect: prefer BETA; else compute from Z + SE; else fall
    # back to OR (b = log(OR)).
    use_beta = beta_i is not None
    use_z = (not use_beta) and (z_i is not None and se_i is not None)
    use_or = (not use_beta) and (not use_z) and (or_i is not None)
    if not (use_beta or use_z or use_or):
        sys.exit("REFUSED: cannot derive signed beta from input (no BETA / "
                 "Z+SE / OR columns).")
    if se_i is None and not use_z:
        sys.exit("REFUSED: cannot derive SE — need SE column or Z column "
                 "(SE inferred from |b|/|Z|).")

    # Pre-load LD's snp.info so we can fill NA freq + sanity-check alleles.
    # SBayesRC::tidy() rejects any row with NA freq (the rate2pq filter
    # divides by NaN), so without this fill every GWAS Catalog harmonised
    # file (which sets hm_effect_allele_frequency=NA) would produce zero
    # SNPs.
    ld_snp_info = _load_ld_snp_info(ld_dir) if ld_dir else {}
    n_freq_filled = 0
    n_allele_mismatch = 0

    out_cojo.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with _open_text(input_path) as f:
        f.readline()  # skip header
        with open(out_cojo, "w") as g:
            g.write("\t".join(COJO_COLS) + "\n")
            for line in f:
                r = line.rstrip("\n").split("\t")
                if len(r) < max(filter(None, (snp_i, a1_i, a2_i, p_i,
                                                beta_i, se_i, z_i, or_i,
                                                frq_i, n_i, n_cas_i,
                                                n_con_i))) + 1:
                    continue
                snp = r[snp_i]; a1 = r[a1_i].upper(); a2 = r[a2_i].upper()
                p = r[p_i]
                if snp in _NA_TOKENS or p in _NA_TOKENS:
                    continue
                # signed effect
                try:
                    if use_beta:
                        b = float(r[beta_i])
                        se_val = float(r[se_i]) if se_i is not None else None
                    elif use_z:
                        z = float(r[z_i])
                        se_val = float(r[se_i])
                        b = z * se_val
                    else:  # OR
                        import math
                        b = math.log(float(r[or_i]))
                        se_val = float(r[se_i])
                except (ValueError, IndexError):
                    continue
                if se_val is None or se_val == 0:
                    continue
                # frequency: prefer input, fall back to LD snp.info A1Freq.
                if frq_i is not None and r[frq_i] not in _NA_TOKENS:
                    frq = r[frq_i]
                else:
                    frq = "NA"
                    ld_entry = ld_snp_info.get(snp)
                    if ld_entry is not None:
                        ld_a1, ld_a2, ld_frq = ld_entry
                        if a1 == ld_a1 and a2 == ld_a2:
                            frq = f"{ld_frq:.6g}"
                            n_freq_filled += 1
                        elif a1 == ld_a2 and a2 == ld_a1:
                            frq = f"{1.0 - ld_frq:.6g}"
                            n_freq_filled += 1
                        else:
                            # Allele set doesn't match; SBayesRC would drop
                            # this anyway. Skip rather than emit NA.
                            n_allele_mismatch += 1
                            continue
                    elif ld_snp_info:
                        # LD has snp.info AND the SNP isn't in HM3 — it'll
                        # be dropped by tidy() at the intersect step. Skip
                        # now to keep the COJO tight.
                        continue
                # N
                if n_i is not None and r[n_i] not in _NA_TOKENS:
                    N_val = r[n_i]
                elif n_cas_i is not None and n_con_i is not None:
                    try:
                        N_val = str(int(r[n_cas_i]) + int(r[n_con_i]))
                    except (ValueError, IndexError):
                        N_val = str(N_const) if N_const else "NA"
                else:
                    N_val = str(N_const) if N_const else "NA"
                if N_val == "NA":
                    continue
                g.write(f"{snp}\t{a1}\t{a2}\t{frq}\t{b:.6g}\t"
                        f"{se_val:.6g}\t{p}\t{N_val}\n")
                n_written += 1
    return n_written


# --------------------------- Apptainer / remote-exec backend ---------------------------

# zhiliz/sbayesrc is x86_64-only. On macOS arm64, two paths:
#   1) Docker/Colima with --platform linux/amd64 (Rosetta emulation; slow)
#   2) SSH to an x86_64 Linux host (e.g. Hoffman2) and run apptainer there
#
# The remote path expects the host to have apptainer/singularity installed
# and the SBayesRC image either already pulled or accessible via
# `docker://zhiliz/sbayesrc`. The skill ships the COJO + LD ref via the
# user's existing remote mount (no extra rsync) so the same paths work on
# both sides.

def _run_sbayesrc_apptainer(cojo_path: Path, ld_dir: Path,
                              annot_cache: Path, out_prefix: str,
                              args) -> int:
    image = "docker://zhiliz/sbayesrc"
    # Build the in-container Rscript invocation. The container ships the R
    # script's logic baked in via `/opt/SBayesRC/run.R`; we mirror our local
    # 3-step pipeline by calling the R API directly.
    r_inline = (
        'suppressPackageStartupMessages(library(SBayesRC)); '
        f'SBayesRC::tidy(mafile="{cojo_path}", LDdir="{ld_dir}", '
        f'output="{out_prefix}.tidy", log2file=TRUE); '
        f'SBayesRC::impute(mafile="{out_prefix}.tidy.ma", '
        f'LDdir="{ld_dir}", output="{out_prefix}.imputed", log2file=TRUE); '
        f'SBayesRC::sbayesrc(mafile="{out_prefix}.imputed.ma", '
        f'LDdir="{ld_dir}", annot="{annot_cache}", '
        f'output="{out_prefix}", log2file=TRUE)'
    )

    if args.remote_host:
        # ssh <host> apptainer exec docker://zhiliz/sbayesrc Rscript -e '...'
        # Assumes paths are reachable on the remote (typically true when the
        # user has the Hoffman2 mount in both places).
        cmd = ["ssh", args.remote_host,
               "apptainer", "exec", image, "Rscript", "-e", r_inline]
        print(f"[prs] backend=apptainer (remote={args.remote_host})  "
              f"{' '.join(cmd[:5])} ...", file=sys.stderr)
        return subprocess.call(cmd)

    # Local execution: prefer apptainer if installed, else Docker with
    # platform emulation.
    if shutil.which("apptainer"):
        cmd = ["apptainer", "exec", image, "Rscript", "-e", r_inline]
    elif shutil.which("docker"):
        # Mount the working dirs the container needs to see.
        mounts = []
        for host_path in {str(cojo_path.parent),
                          str(ld_dir),
                          str(annot_cache.parent),
                          str(Path(out_prefix).parent)}:
            mounts += ["-v", f"{host_path}:{host_path}"]
        cmd = ["docker", "run", "--rm", "--platform", "linux/amd64",
               *mounts, "zhiliz/sbayesrc",
               "Rscript", "-e", r_inline]
        print(f"[prs] backend=apptainer (via Docker w/ amd64 emulation)",
              file=sys.stderr)
    else:
        sys.exit("REFUSED: --backend apptainer needs either `apptainer` or "
                 "`docker` on PATH (or --remote-host <ssh-target>). "
                 "Install Apptainer/Singularity, install Docker, or run on "
                 "a Linux x86_64 host.")
    return subprocess.call(cmd)


# --------------------------- Method dispatch ---------------------------

def run_sbayesrc(args, panel: dict, ld_dir: Path, cojo_path: Path,
                  okg_node_ids: dict) -> int:
    backend = args.backend
    out_prefix = str(args.out)
    annot_cache = CACHE_ROOT / "annot_baseline2.2.txt"
    if backend == "r":
        r_script = Path(__file__).parent / "sbayesrc_run.R"
        if not r_script.exists():
            sys.exit(f"ERROR: companion R script not found at {r_script}; "
                     f"the prs skill is mis-installed.")
        cmd = ["Rscript", str(r_script),
               "--cojo", str(cojo_path),
               "--ld-dir", str(ld_dir),
               "--annot-cache", str(annot_cache),
               "--out", out_prefix]
        print(f"[prs] backend=r  {' '.join(cmd)}", file=sys.stderr)
        rc = subprocess.call(cmd)
    elif backend == "apptainer":
        rc = _run_sbayesrc_apptainer(cojo_path, ld_dir, annot_cache,
                                       out_prefix, args)
    else:
        sys.exit(f"REFUSED: unknown backend {backend!r}; "
                 f"allowed: r, apptainer")

    # SBayesRC writes per-SNP weights to <prefix>.txt and posterior fit
    # statistics to <prefix>.par. Earlier versions used .snpRes; we accept
    # either.
    weights_path = Path(f"{out_prefix}.txt")
    if not weights_path.exists():
        legacy = Path(f"{out_prefix}.snpRes")
        if legacy.exists():
            weights_path = legacy
    fit = _parse_sbayesrc_par(Path(f"{out_prefix}.par"))
    summary = {"weights_path": str(weights_path)}
    if weights_path.exists():
        try:
            with open(weights_path) as f:
                summary["n_snps_retained"] = sum(1 for _ in f) - 1
        except Exception:
            pass

    manifest = {
        "method": "sbayesrc",
        "gwas_sumstats_input": str(args.gwas_sumstats),
        "output_prefix": out_prefix,
        "weights_path": str(weights_path),
        "n_snps_retained": summary.get("n_snps_retained"),
        "okg_node_ids": okg_node_ids,
        "ld_reference": {
            "node_id": panel["node_id"],
            "local_path": str(ld_dir),
            "source_url": panel.get("source_url"),
            "genome_build": panel.get("genome_build"),
            "ancestry": panel.get("ancestry"),
        },
        "method_specific": fit,
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    mpath = Path(f"{out_prefix}.prs.json")
    mpath.parent.mkdir(parents=True, exist_ok=True)
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"manifest -> {mpath}", file=sys.stderr)
    if rc == 0:
        hsq = fit.get("hsq")
        hsq_se = fit.get("hsq_se")
        hsq_str = f"hsq={hsq:.3f} (SE {hsq_se:.3f})" if hsq is not None else ""
        print(f"sbayesrc weights -> {weights_path} "
              f"({summary.get('n_snps_retained', '?')} SNPs)  {hsq_str}")
    return rc


def _parse_sbayesrc_par(par: Path) -> dict:
    """Parse SBayesRC's `<prefix>.par` posterior summary file.

    Format is `Item<TAB>Mean<TAB>SD` with rows like:
        hsq    0.2586    0.0176
        nnz    4042.5    705.7
        SigmaSq 0.00124  0.000113
        ...
    Returns a dict keyed on lowercased item name with `_se` siblings for SD.
    """
    out = {}
    if not par.exists():
        return out
    try:
        with open(par) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                name = parts[0].lower()
                try:
                    out[name] = float(parts[1])
                    out[f"{name}_se"] = float(parts[2])
                except ValueError:
                    continue
    except Exception:
        pass
    return out


def _parse_sbayesrc_log(log: Path) -> dict:
    text = log.read_text(errors="ignore")
    out = {}
    for key, pat in [
        ("hsq", r"hsq\s*=\s*([-\d.eE+]+)"),
        ("hsq_se", r"hsq[_\s]se\s*=\s*([-\d.eE+]+)"),
        ("polygenicity_pi", r"Pi\s*=\s*([-\d.eE+]+)"),
        ("n_mcmc_iter", r"niter\s*=\s*(\d+)"),
        ("n_burnin", r"nburn\s*=\s*(\d+)"),
    ]:
        m = re.search(pat, text)
        if m:
            out[key] = float(m.group(1)) if "." in m.group(1) or "e" in m.group(1).lower() else int(m.group(1))
    return out


# --------------------------- CLI ---------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--method", type=str, default="sbayesrc",
                   choices=["sbayesrc"],
                   help="PRS method to run (default: sbayesrc). More "
                        "methods can be wired in by adding a branch below.")
    p.add_argument("--backend", type=str, default="r",
                   choices=["r", "apptainer"],
                   help="Execution backend for sbayesrc. `r` uses the local "
                        "R package install (default, works on macOS arm64 "
                        "with the Makevars patch). `apptainer` runs the "
                        "x86_64-only docker://zhiliz/sbayesrc image via "
                        "apptainer/singularity locally or via Docker with "
                        "linux/amd64 emulation; use --remote-host to "
                        "route to a Linux cluster (e.g. Hoffman2).")
    p.add_argument("--remote-host", dest="remote_host", type=str,
                   default=None,
                   help="When --backend apptainer is set, ssh target where "
                        "apptainer is installed. Paths in COJO + LD dir + "
                        "out prefix must be visible on the remote (usually "
                        "via a shared mount).")
    p.add_argument("--gwas-sumstats", dest="gwas_sumstats", type=Path,
                   required=True,
                   help="Path to GWAS sumstats (raw harmonised TSV, LDSC "
                        ".sumstats.gz, or pre-converted COJO TSV).")
    p.add_argument("--out", type=Path, required=True,
                   help="Output prefix.")
    p.add_argument("--ancestry", type=str,
                   choices=["eur", "eas", "afr"],
                   help="Ancestry of the GWAS. Picks the matching SBayesRC "
                        "HM3 LD reference. Required unless "
                        "--okg-ld-panel-id is given.")
    p.add_argument("--okg-ld-panel-id", dest="okg_ld_panel_id", type=str,
                   help="Explicit OKG ld_panel node ID (overrides --ancestry).")
    p.add_argument("--okg-dataset-id", dest="okg_dataset_id", type=str,
                   help="OKG dataset_metadata node to auto-resolve N.")
    p.add_argument("--N", type=int, default=None,
                   help="Constant total sample size (overrides OKG-resolved).")
    p.add_argument("--repo-cache", type=Path, default=CACHE_ROOT,
                   help=f"Cache root for LD references (default: {CACHE_ROOT})")
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None,
                   help="OKG repo for metadata resolution (honors $OKG_REPO).")
    args = p.parse_args()

    panel = resolve_ld_panel(args.okg_repo, args.ancestry,
                              args.okg_ld_panel_id)
    N_resolved = args.N or resolve_dataset_n(args.okg_repo,
                                              args.okg_dataset_id)

    okg_node_ids = dict(OKG_NODES_SBAYESRC)
    okg_node_ids["ld_panel"] = panel["node_id"]
    if args.okg_dataset_id:
        okg_node_ids["dataset"] = args.okg_dataset_id

    ld_dir = ensure_ld_dir(panel, Path(args.repo_cache))
    print(f"[prs] LD reference: {ld_dir} (panel={panel['node_id']}, "
          f"build={panel.get('genome_build')}, "
          f"ancestry={panel.get('ancestry')})", file=sys.stderr)

    cojo_path = Path(str(args.out) + ".cojo.tsv")
    n_written = detect_format_and_convert_to_cojo(
        args.gwas_sumstats, N_resolved, cojo_path, ld_dir=ld_dir)
    print(f"[prs] converted to COJO: {cojo_path} ({n_written:,} rows)",
          file=sys.stderr)
    if n_written < 50_000:
        sys.exit(f"REFUSED: COJO has {n_written:,} rows (<50,000). "
                 f"SBayesRC needs the full HapMap3 scaffold.")

    if args.method == "sbayesrc":
        return run_sbayesrc(args, panel, ld_dir, cojo_path, okg_node_ids)
    # Unreachable given `choices=["sbayesrc"]`, but the dispatch is here for
    # future methods.
    sys.exit(f"REFUSED: unknown method {args.method!r}")


if __name__ == "__main__":
    sys.exit(main())
