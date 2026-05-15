#!/usr/bin/env python3
"""twas.py --- orchestrator for the `twas` skill.

Runs TWAS-FUSION (Gusev et al. 2016) against pre-computed GTEx v8
multi-tissue weight panels. Resolves the panel via the OKG when
`$OKG_REPO` is set; downloads + caches the panel and FUSION 1000G LDREF
on first use; normalizes sumstats to FUSION's expected
SNP/A1/A2/Z layout; runs FUSION.assoc_test.R per chromosome; writes a
sidecar manifest citing the method/software/paper/ld_panel/dataset/
tissue/cohort OKG nodes.

Usage:
    twas.py --sumstats <file> --tissue Brain_Cortex --ancestry EUR --out twas/AD_BrainCortex
    twas.py --sumstats <file> --okg-dataset-id dataset:fusion_gtex_v8_eur:liver --out plots/T01_liver
    twas.py --sumstats <file> --okg-trait-id trait:height --out twas/height_multi
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

CACHE_ROOT = Path(os.environ.get("TWAS_CACHE",
                                  str(Path.home() / ".cache" / "twas")))
FUSION_REPO_URL = "https://github.com/gusevlab/fusion_twas.git"
FUSION_LDREF_URL = "https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2"
HERE = Path(__file__).resolve().parent

EXIT_BAD_INPUT = 2
EXIT_PANEL_UNRESOLVED = 3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------- sumstats normalisation -----------------------

def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def detect_columns(header: list[str]) -> dict[str, int]:
    """Map a sumstats header to FUSION-needed column indices."""
    cols = {c: i for i, c in enumerate(header)}
    cols_l = {c.lower(): i for i, c in enumerate(header)}

    def first(*candidates):
        for c in candidates:
            if c in cols:
                return cols[c]
            if c.lower() in cols_l:
                return cols_l[c.lower()]
        return None

    return {
        "snp":  first("SNP", "snp", "hm_rsid", "rsid", "MarkerName", "markername"),
        "a1":   first("A1", "a1", "Effect_allele", "effect_allele", "hm_effect_allele", "Allele1"),
        "a2":   first("A2", "a2", "Other_allele", "other_allele",
                       "Non_Effect_allele", "hm_other_allele", "Allele2"),
        "z":    first("Z", "z", "zscore"),
        "beta": first("BETA", "beta", "Beta", "hm_beta", "Effect", "b"),
        "se":   first("SE", "se", "standard_error", "StdErr", "stderr"),
        "or":   first("OR", "or", "odds_ratio", "hm_odds_ratio"),
        "p":    first("P", "p", "Pvalue", "p_value", "pval"),
    }


def normalise_to_fusion(input_path: Path, work_dir: Path) -> tuple[Path, dict]:
    """Read sumstats, write FUSION-style 'SNP A1 A2 Z' (whitespace-delim) TSV.

    Returns (out_path, info_dict with n_rows / n_kept / n_with_z / sep).
    """
    out = work_dir / "fusion_sumstats.tsv"
    sep_used = "\t"
    n_kept = n_rows = 0
    with _open_text(input_path) as fin:
        first_line = fin.readline().rstrip("\n")
        sep = "\t" if "\t" in first_line else " "
        sep_used = sep
        header = first_line.split(sep)
        idx = detect_columns(header)
        if idx["snp"] is None or idx["a1"] is None or idx["a2"] is None:
            sys.exit(f"REFUSED: sumstats missing SNP/A1/A2 (looked for "
                     f"SNP/hm_rsid/MarkerName, A1/Effect_allele/Allele1, "
                     f"A2/Other_allele/Allele2/Non_Effect_allele); header = {header[:25]}")
        have_z = idx["z"] is not None
        have_betase = idx["beta"] is not None and idx["se"] is not None
        have_orse = idx["or"] is not None and idx["se"] is not None
        if not (have_z or have_betase or have_orse):
            sys.exit("REFUSED: sumstats has neither Z nor BETA+SE nor OR+SE; "
                     "can't construct Z.")
        with open(out, "w") as fout:
            fout.write("SNP\tA1\tA2\tZ\n")
            for line in fin:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split(sep)
                if len(fields) < len(header):
                    continue
                n_rows += 1
                snp = fields[idx["snp"]] if idx["snp"] is not None else None
                a1 = fields[idx["a1"]] if idx["a1"] is not None else None
                a2 = fields[idx["a2"]] if idx["a2"] is not None else None
                if not snp or snp in ("NA", ".") or not a1 or not a2:
                    continue
                z = None
                if have_z:
                    try:
                        z = float(fields[idx["z"]])
                    except (ValueError, IndexError):
                        z = None
                if z is None and have_betase:
                    try:
                        b = float(fields[idx["beta"]])
                        s = float(fields[idx["se"]])
                        if s > 0:
                            z = b / s
                    except (ValueError, IndexError):
                        z = None
                if z is None and have_orse:
                    try:
                        o = float(fields[idx["or"]])
                        s = float(fields[idx["se"]])
                        if o > 0 and s > 0:
                            z = math.log(o) / s
                    except (ValueError, IndexError):
                        z = None
                if z is None or not math.isfinite(z):
                    continue
                # rsids are conventionally lowercase ("rs8100066"); only
                # the alleles get uppercased to match panel/LDREF style.
                fout.write(f"{snp}\t{a1.upper()}\t{a2.upper()}\t{z:.6g}\n")
                n_kept += 1
    info = {"n_rows": n_rows, "n_kept": n_kept, "sep_in": sep_used}
    sys.stderr.write(f"normalised: kept {n_kept:,} of {n_rows:,} rows -> {out}\n")
    return out, info


# ---------------------------- OKG lookup -----------------------------------

def okg_mcp_call(okg_repo: Path, tool: str, arguments: dict) -> dict | None:
    """Call a tool on the statgen-analysis MCP server via stdio."""
    server = okg_repo / "deployments/statgen-analysis/server.py"
    if not server.exists():
        return None
    env = os.environ.copy()
    env.setdefault("OKG_DSN",
                   "postgres://postgres:okg@localhost:5449/statgen_analysis")
    try:
        proc = subprocess.Popen(
            ["uv", "run", "--extra", "mcp", "python", str(server)],
            cwd=str(okg_repo), env=env,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True, bufsize=1,
        )
    except FileNotFoundError:
        return None
    def send(m): proc.stdin.write(json.dumps(m) + "\n"); proc.stdin.flush()
    def read(): return json.loads(proc.stdout.readline())
    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                         "clientInfo": {"name": "twas", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": tool, "arguments": arguments}})
        resp = read()
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    sc = resp.get("result", {}).get("structuredContent") or {}
    return sc


def tissue_slug(panel: str) -> str:
    s = panel.lower().replace("-", "_")
    return re.sub(r"_+", "_", s).strip("_")


def resolve_panel(args) -> tuple[str, dict]:
    """Return (dataset_node_id, attrs). Either via OKG MCP or by URL-pattern fallback."""
    if args.okg_dataset_id:
        if args.okg_repo:
            attrs = okg_mcp_call(args.okg_repo, "get_node",
                                  {"node_id": args.okg_dataset_id})
            if attrs and attrs.get("attrs"):
                return args.okg_dataset_id, attrs["attrs"]
        # Fallback: derive URL from node id pattern
        m = re.match(r"^dataset:fusion_gtex_v8_(eur|all):(.+)$", args.okg_dataset_id)
        if not m:
            sys.exit(f"REFUSED: --okg-dataset-id doesn't look like a FUSION "
                     f"GTEx v8 panel: {args.okg_dataset_id}")
        anc, slug = m.group(1).upper(), m.group(2)
        # We need the FUSION tissue token (e.g. Brain_Cortex), which we can't
        # exactly reverse from the slug — require --tissue alongside if no OKG.
        if not args.tissue:
            sys.exit("REFUSED: without OKG_REPO, --okg-dataset-id requires "
                     "--tissue to reconstruct the FUSION panel name.")
        url = (f"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/{anc}/"
               f"GTExv8.{anc}.{args.tissue}.tar.gz")
        return args.okg_dataset_id, {
            "dataset_id": args.okg_dataset_id, "source_url": url,
            "tissue": args.tissue, "ancestry_scope": "european" if anc == "EUR" else "multi_ancestry",
        }
    if args.tissue and args.ancestry:
        ancestry_tag = args.ancestry.lower()
        slug = tissue_slug(args.tissue)
        node_id = f"dataset:fusion_gtex_v8_{ancestry_tag}:{slug}"
        url = (f"https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/{args.ancestry}/"
               f"GTExv8.{args.ancestry}.{args.tissue}.tar.gz")
        attrs = {"dataset_id": node_id, "source_url": url,
                 "tissue": args.tissue,
                 "ancestry_scope": "european" if args.ancestry == "EUR" else "multi_ancestry"}
        if args.okg_repo:
            okg_attrs = okg_mcp_call(args.okg_repo, "get_node",
                                      {"node_id": node_id})
            if okg_attrs and okg_attrs.get("attrs"):
                attrs = okg_attrs["attrs"]
        return node_id, attrs
    sys.exit("REFUSED: provide --tissue + --ancestry, or --okg-dataset-id")


# ---------------------------- downloads + caches ---------------------------

def ensure_fusion_repo() -> Path:
    repo_dir = CACHE_ROOT / "repo"
    if (repo_dir / "FUSION.assoc_test.R").exists():
        return repo_dir
    repo_dir.mkdir(parents=True, exist_ok=True)
    sys.stderr.write(f"cloning FUSION repo -> {repo_dir}\n")
    subprocess.check_call(["git", "clone", "--depth=1", FUSION_REPO_URL, str(repo_dir)])
    return repo_dir


def ensure_ldref() -> Path:
    ldref_dir = CACHE_ROOT / "ldref"
    if (ldref_dir / "LDREF" / "1000G.EUR.1.bed").exists():
        return ldref_dir / "LDREF"
    ldref_dir.mkdir(parents=True, exist_ok=True)
    tarball = ldref_dir / "LDREF.tar.bz2"
    if not tarball.exists():
        sys.stderr.write(f"downloading FUSION LDREF -> {tarball}\n")
        urllib.request.urlretrieve(FUSION_LDREF_URL, tarball)
    sys.stderr.write(f"unpacking LDREF...\n")
    with tarfile.open(tarball, "r:bz2") as tar:
        tar.extractall(ldref_dir)
    return ldref_dir / "LDREF"


def ensure_panel(node_id: str, attrs: dict) -> Path:
    """Download + unpack the FUSION panel; return the unpacked WEIGHTS dir."""
    # Tag from node id like dataset:fusion_gtex_v8_eur:brain_cortex
    panel_tag = node_id.replace(":", "_").replace("/", "_")
    panel_dir = CACHE_ROOT / "panels" / panel_tag
    if any(panel_dir.glob("*.pos")):
        return panel_dir
    url = attrs.get("source_url")
    if not url:
        sys.exit(f"REFUSED: panel node {node_id} has no source_url attr")
    panel_dir.mkdir(parents=True, exist_ok=True)
    tarball = panel_dir / Path(url).name
    if not tarball.exists():
        sys.stderr.write(f"downloading panel -> {tarball}\n")
        urllib.request.urlretrieve(url, tarball)
    sys.stderr.write(f"unpacking panel ({tarball.stat().st_size/1e6:.1f} MB)...\n")
    if tarball.suffix in (".gz",) or tarball.name.endswith(".tar.gz"):
        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(panel_dir)
    elif tarball.name.endswith(".tar.bz2"):
        with tarfile.open(tarball, "r:bz2") as tar:
            tar.extractall(panel_dir)
    else:
        sys.exit(f"unknown archive type: {tarball.name}")
    return panel_dir


def find_pos_file(panel_dir: Path, use_filtered: bool) -> Path:
    """Locate the .pos manifest the FUSION panel ships. Filtered = sig-heritability genes.

    GTEx v8 archives ship `<panel>.pos` (filtered, sig-heritability) and
    `<panel>.nofilter.pos` (all genes). Older FUSION panels use `no_filter`
    (with underscore); accept both forms.
    """
    pos_files = list(panel_dir.rglob("*.pos"))
    if not pos_files:
        sys.exit(f"no .pos files in panel dir {panel_dir}")
    def is_unfiltered(p: Path) -> bool:
        n = p.name.lower()
        return "nofilter" in n or "no_filter" in n
    if use_filtered:
        pref = [p for p in pos_files if not is_unfiltered(p)]
        if pref:
            return pref[0]
    else:
        pref = [p for p in pos_files if is_unfiltered(p)]
        if pref:
            return pref[0]
    return pos_files[0]


# ---------------------------- run FUSION per chromosome ---------------------

def chromosomes_in_pos(pos_file: Path) -> list[int]:
    """Extract the set of chromosomes that have weights in this panel."""
    chrs: set[int] = set()
    with open(pos_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            chr_idx = header.index("CHR")
        except ValueError:
            return list(range(1, 23))
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > chr_idx:
                try:
                    chrs.add(int(parts[chr_idx]))
                except ValueError:
                    continue
    return sorted(chrs)


def run_fusion(repo_dir: Path, sumstats: Path, pos_file: Path,
                panel_dir: Path, ldref_dir: Path, chrom: int,
                out_dat: Path, coloc_p: float | None, perm: int | None,
                log: Path) -> bool:
    # FUSION.assoc_test.R uses here::here("utils","plink_utils.R") so it MUST
    # run with cwd == the fusion_twas repo root (where utils/ lives and
    # where here::here() will anchor). We pass absolute paths for all other
    # files so the run is location-independent except for the cwd anchor.
    assoc_r = "FUSION.assoc_test.R"
    cmd = ["Rscript", assoc_r,
           "--sumstats", str(sumstats.resolve()),
           "--weights", str(pos_file.resolve()),
           "--weights_dir", str(panel_dir.resolve()),
           "--ref_ld_chr", str((ldref_dir / "1000G.EUR.").resolve()),
           "--chr", str(chrom),
           "--out", str(out_dat.resolve())]
    if coloc_p is not None:
        cmd += ["--coloc_P", str(coloc_p)]
    if perm is not None:
        cmd += ["--perm", str(perm)]
    sys.stderr.write(f"FUSION chr{chrom}: cwd={repo_dir}; {' '.join(cmd[:4])}...\n")
    with open(log, "ab") as logf:
        logf.write(f"\n## chr{chrom} ##\n".encode())
        logf.flush()
        rc = subprocess.call(cmd, stdout=logf, stderr=logf,
                              cwd=str(repo_dir))
    return rc == 0 and out_dat.exists()


# ---------------------------- main -----------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sumstats", type=Path, required=True)
    p.add_argument("--tissue", type=str,
                   help="FUSION tissue panel name (e.g. Brain_Cortex)")
    p.add_argument("--ancestry", type=str, choices=["EUR", "ALL"],
                   help="Panel ancestry: EUR-only or ALL-ancestry. Default EUR if --tissue given.")
    p.add_argument("--okg-dataset-id", type=str,
                   help="Explicit OKG dataset_metadata node id for the panel.")
    p.add_argument("--okg-trait-id", type=str,
                   help="OKG trait node id; if set without --tissue, runs over all "
                        "tissues with `relevant_to <trait>` edges.")
    p.add_argument("--okg-gwas-dataset-id", type=str,
                   help="OKG dataset node for the GWAS (recorded in sidecar).")
    p.add_argument("--no-filter", action="store_true",
                   help="Use the panel's no_filter.pos (all genes) instead of "
                        "the filtered sig-heritability .pos (default).")
    p.add_argument("--chr", type=str,
                   help="Comma-separated chromosomes (default: all in panel).")
    p.add_argument("--coloc-p", type=float, default=None)
    p.add_argument("--perm", type=int, default=None)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None)
    args = p.parse_args()

    if args.tissue and not args.ancestry:
        args.ancestry = "EUR"

    args.sumstats = args.sumstats.expanduser().resolve()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if not args.sumstats.exists():
        sys.exit(f"sumstats not found: {args.sumstats}")

    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix="twas_"))

    # --- Set up dependencies (cached after first use) ---
    repo_dir = ensure_fusion_repo()
    ldref_dir = ensure_ldref()

    # --- Resolve the panel ---
    node_id, attrs = resolve_panel(args)
    sys.stderr.write(f"panel: {node_id}  url={attrs.get('source_url')}\n")
    panel_dir = ensure_panel(node_id, attrs)
    pos_file = find_pos_file(panel_dir, use_filtered=not args.no_filter)
    sys.stderr.write(f"pos: {pos_file}\n")

    # --- Normalise sumstats ---
    fusion_ss, ss_info = normalise_to_fusion(args.sumstats, work_dir)

    # --- Determine chromosomes to run ---
    if args.chr:
        chroms = [int(x) for x in args.chr.split(",")]
    else:
        chroms = chromosomes_in_pos(pos_file)
    sys.stderr.write(f"running FUSION on {len(chroms)} chromosomes: {chroms}\n")

    # --- Run per chromosome ---
    log = Path(f"{args.out}.log")
    log.write_text(f"twas run at {dt.datetime.now(dt.timezone.utc).isoformat()}\n"
                    f"panel: {node_id}\nsumstats: {args.sumstats}\n")
    out_assoc = Path(f"{args.out}.assoc.tsv")
    per_chr_paths: list[Path] = []
    failures: list[int] = []
    for chrom in chroms:
        per = work_dir / f"chr{chrom}.dat"
        ok = run_fusion(repo_dir, fusion_ss, pos_file, panel_dir, ldref_dir,
                         chrom, per, args.coloc_p, args.perm, log)
        if ok:
            per_chr_paths.append(per)
        else:
            failures.append(chrom)
    if not per_chr_paths:
        sys.exit(f"FUSION produced no output for any chromosome; see {log}")

    # --- Concatenate ---
    header_written = False
    with open(out_assoc, "w") as fout:
        for p_ in per_chr_paths:
            with open(p_) as fin:
                hdr = fin.readline()
                if not header_written:
                    fout.write(hdr)
                    header_written = True
                shutil.copyfileobj(fin, fout)

    # --- Parse + summarise ---
    n_genes = 0
    sig = []
    top = []
    with open(out_assoc) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            id_i = header.index("ID")
            z_i = header.index("TWAS.Z")
            p_i = header.index("TWAS.P")
            chr_i = header.index("CHR")
            eqtl_i = header.index("BEST.GWAS.ID") if "BEST.GWAS.ID" in header else None
        except ValueError as e:
            sys.exit(f"FUSION output missing expected column: {e}; header={header}")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(id_i, z_i, p_i, chr_i):
                continue
            try:
                tp = float(fields[p_i]) if fields[p_i] != "NA" else None
                tz = float(fields[z_i]) if fields[z_i] != "NA" else None
            except ValueError:
                continue
            if tp is None or tz is None:
                continue
            n_genes += 1
            entry = {
                "id": fields[id_i], "chr": int(fields[chr_i]),
                "twas_z": tz, "twas_p": tp,
                "best_eqtl": fields[eqtl_i] if eqtl_i is not None else None,
            }
            top.append(entry)
    bonf = 0.05 / n_genes if n_genes > 0 else 1
    sig = [t for t in top if t["twas_p"] < bonf]
    top_sorted = sorted(top, key=lambda x: x["twas_p"])[:10]

    # --- Manifest ---
    okg_node_ids = {
        "method": "method:fusion_twas",
        "software": "software:fusion",
        "paper": "paper:fusion_2016",
        "ld_panel": "ld_panel:fusion_1000g_eur",
        "cohort": "cohort:gtex_v8",
        "dataset_panel": node_id,
    }
    if args.tissue:
        okg_node_ids["tissue"] = f"tissue:{tissue_slug(args.tissue)}"
    if args.okg_gwas_dataset_id:
        okg_node_ids["dataset_gwas"] = args.okg_gwas_dataset_id

    # FUSION repo commit
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        commit = None

    manifest = {
        "sumstats_input": str(args.sumstats),
        "sumstats_sha256": sha256(args.sumstats),
        "sumstats_normalise": ss_info,
        "output_assoc": str(out_assoc),
        "output_log": str(log),
        "tissue": args.tissue,
        "ancestry": args.ancestry,
        "panel_url": attrs.get("source_url"),
        "fusion_repo": FUSION_REPO_URL,
        "fusion_commit": commit,
        "ldref_url": FUSION_LDREF_URL,
        "use_filtered_pos": not args.no_filter,
        "pos_file": str(pos_file),
        "chromosomes_run": chroms,
        "chromosomes_failed": failures,
        "n_genes_tested": n_genes,
        "n_significant_bonferroni": len(sig),
        "bonferroni_threshold": bonf,
        "top_hits": top_sorted,
        "okg_node_ids": okg_node_ids,
        "captured_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    sidecar = Path(f"{args.out}.twas.json")
    sidecar.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"twas results -> {out_assoc}")
    print(f"  n_genes_tested = {n_genes}")
    print(f"  n_significant (p < {bonf:.2e}) = {len(sig)}")
    print(f"  top hit: {top_sorted[0]['id']} (TWAS.P={top_sorted[0]['twas_p']:.2e}, "
          f"TWAS.Z={top_sorted[0]['twas_z']:.2f})" if top_sorted else "")
    print(f"sidecar -> {sidecar}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
