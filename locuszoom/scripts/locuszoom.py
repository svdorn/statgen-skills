#!/usr/bin/env python3
"""locuszoom.py --- orchestrator for the `locuszoom` skill.

Renders a LocusZoom-style regional plot from a GWAS sumstats file using
the locuszoomr R package + LDlink REST API.

Two-tier token resolution: --ldlink-token (and cache), then $LDLINK_TOKEN,
then ~/.cache/locuszoom/ldlink_token. Exits with `LDLINK_TOKEN_MISSING`
if none found so the agent can prompt the user.

Usage:
    locuszoom.py --sumstats <file> --gene APOE --out plots/AD_APOE
    locuszoom.py --sumstats <file> --region 19:45000000-45500000 --out plots/X
    locuszoom.py --sumstats <file> --lead-snp rs429358 --flank 200000 --out plots/Y
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

CACHE_ROOT = Path(os.environ.get("LOCUSZOOM_CACHE",
                                  str(Path.home() / ".cache" / "locuszoom")))
TOKEN_FILE = CACHE_ROOT / "ldlink_token"
HERE = Path(__file__).resolve().parent

EXIT_TOKEN_MISSING = 2
EXIT_TOKEN_INVALID = 3
EXIT_BAD_LOCUS = 4
EXIT_SPARSE_LOCUS = 5


def resolve_token(cli_token: str | None) -> tuple[str, str]:
    """Return (token, source). Caches cli_token if supplied."""
    if cli_token:
        CACHE_ROOT.mkdir(parents=True, exist_ok=True)
        TOKEN_FILE.write_text(cli_token.strip() + "\n")
        try:
            os.chmod(TOKEN_FILE, 0o600)
        except OSError:
            pass
        return cli_token.strip(), "cli_flag_cached"
    env = os.environ.get("LDLINK_TOKEN")
    if env:
        return env.strip(), "env_LDLINK_TOKEN"
    if TOKEN_FILE.exists():
        token = TOKEN_FILE.read_text().strip()
        if token:
            return token, f"cache:{TOKEN_FILE}"
    sys.stderr.write(
        "LDLINK_TOKEN_MISSING\n"
        "No LDlink API token found. Get one (free, ~30 seconds) at:\n"
        "  https://ldlink.nih.gov/?tab=apiaccess\n"
        f"Then re-invoke with `--ldlink-token <token>` (it caches to {TOKEN_FILE}).\n"
    )
    sys.exit(EXIT_TOKEN_MISSING)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def detect_columns(header: list[str]) -> dict[str, str]:
    """Map a sumstats header to standard column names.

    Returns dict with keys rsid, chrom, pos, p, beta, se. Values are the
    actual header strings. Raises if required ones are missing.
    """
    names = {h.lower(): h for h in header}
    def first(*candidates):
        for c in candidates:
            if c in names:
                return names[c]
        return None
    out = {
        "rsid":  first("hm_rsid", "rsid", "snp", "markername", "rs_id"),
        "chrom": first("hm_chrom", "chromosome", "chr", "#chrom", "#chr"),
        "pos":   first("hm_pos", "base_pair_location", "bp", "position", "pos"),
        "p":     first("p_value", "p", "pvalue", "p_val"),
        "beta":  first("hm_beta", "beta", "effect", "b"),
        "se":    first("standard_error", "se", "stderr"),
    }
    missing = [k for k in ("rsid", "p") if out[k] is None]
    if missing:
        raise SystemExit(f"sumstats missing required columns: {missing}; "
                         f"found header = {header[:20]}")
    return out


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def normalise_sumstats(input_path: Path, work_dir: Path) -> Path:
    """Read input sumstats, write a normalised TSV (rsid, chrom, pos, p, beta, se).

    Output is unzipped TSV at <work_dir>/normalised.tsv.
    """
    out = work_dir / "normalised.tsv"
    with _open(input_path) as fin:
        first = fin.readline().rstrip("\n")
        sep = "\t" if "\t" in first else (" " if " " in first else "\t")
        header = first.split(sep)
        cols = detect_columns(header)
        idx = {k: header.index(v) for k, v in cols.items() if v is not None}
        with open(out, "w") as fout:
            fout.write("rsid\tchrom\tpos\tp\tbeta\tse\n")
            n_rows = 0
            n_kept = 0
            for line in fin:
                fields = line.rstrip("\n").split(sep)
                if len(fields) < len(header):
                    continue
                n_rows += 1
                rsid = fields[idx["rsid"]] if "rsid" in idx else "NA"
                ch = fields[idx["chrom"]] if "chrom" in idx else "NA"
                ps = fields[idx["pos"]] if "pos" in idx else "NA"
                p  = fields[idx["p"]] if "p" in idx else "NA"
                b  = fields[idx["beta"]] if "beta" in idx else "NA"
                s  = fields[idx["se"]] if "se" in idx else "NA"
                if p in ("", "NA", "nan") or rsid in ("", "NA", "."):
                    continue
                n_kept += 1
                fout.write(f"{rsid}\t{ch}\t{ps}\t{p}\t{b}\t{s}\n")
        sys.stderr.write(f"normalised: kept {n_kept:,} of {n_rows:,} rows -> {out}\n")
    return out


def parse_region(region: str) -> tuple[str, int, int]:
    m = re.match(r"^([0-9XYM]+|chr[0-9XYM]+):(\d+)[-:](\d+)$", region)
    if not m:
        sys.stderr.write(f"BAD_REGION: --region must be like '19:45000000-45500000', got {region!r}\n")
        sys.exit(EXIT_BAD_LOCUS)
    chrom = m.group(1).replace("chr", "")
    start, end = int(m.group(2)), int(m.group(3))
    if start >= end:
        sys.stderr.write(f"BAD_REGION: start ({start}) must be < end ({end})\n")
        sys.exit(EXIT_BAD_LOCUS)
    return chrom, start, end


def _materialize_rsid_pip(cs_path: Path, gwas_extras_path: Path,
                           work_dir: Path) -> Path:
    """The finemap skill's `<prefix>.sushie.cs.tsv` keys variants by `chr:pos`,
    not rsid; build an rsid-keyed PIP TSV by joining against
    `<prefix>.gwas_extras.tsv` (chrpos → rsid map). Returns the joined path."""
    rsid_map: dict[str, str] = {}
    with open(gwas_extras_path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                rsid_map[parts[0]] = parts[1]
    out = work_dir / "finemap_cs_rsid.tsv"
    n_in = n_out = 0
    with open(cs_path) as fin, open(out, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        try:
            cs_i = header.index("CSIndex")
            snp_i = header.index("snp")
            pip_i = header.index("pip_all")
        except ValueError:
            sys.exit(f"finemap cs file {cs_path} missing CSIndex/snp/pip_all "
                     f"columns; header = {header}")
        fout.write("snp\tpip\tCSIndex\n")
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            n_in += 1
            chrpos = fields[snp_i]
            rsid = rsid_map.get(chrpos)
            if rsid is None:
                continue
            n_out += 1
            fout.write(f"{rsid}\t{fields[pip_i]}\t{fields[cs_i]}\n")
    sys.stderr.write(f"finemap rsid-mapped: {n_out}/{n_in} CS variants -> {out}\n")
    return out


def resolve_finemap_pip(args, work_dir: Path) -> tuple[Path | None, dict]:
    """Return (pip_tsv_path, finemap_metadata). Either `--finemap-pip` directly,
    or `--finemap-sidecar` → derived from a `finemap` skill `.finemap.json`.
    Returns (None, {}) if neither flag was given.

    For `--finemap-sidecar`, supports three modes (in priority order):
      1. Sidecar carries an explicit `cs_path` / `weights_path` / `weight_path` /
         `pip_file` key that points at an existing file.
      2. Sidecar carries `output_prefix` (the sushie finemap skill convention):
         look for `<output_prefix>.sushie.cs.tsv` and join against
         `<output_prefix>.gwas_extras.tsv` to translate chr:pos → rsid.
      3. Sidecar carries `output_prefix` and `<output_prefix>.sushie.weights.tsv`
         exists: use its `sushie_pip_all` column directly (no rsid join — works
         if upstream sumstats were already rsid-keyed)."""
    meta: dict = {}
    if args.finemap_pip:
        path = args.finemap_pip.expanduser().resolve()
        if not path.exists():
            sys.exit(f"--finemap-pip not found: {path}")
        meta = {"source": "finemap_pip_direct"}
        return path, meta
    if args.finemap_sidecar:
        sc_path = args.finemap_sidecar.expanduser().resolve()
        if not sc_path.exists():
            sys.exit(f"--finemap-sidecar not found: {sc_path}")
        sc = json.loads(sc_path.read_text())
        okg_passthrough = {
            "okg_method": sc.get("okg_node_ids", {}).get("method"),
            "okg_software": (sc.get("okg_node_ids", {}).get("software")
                              or sc.get("okg_node_ids", {}).get("software_canonical")
                              or sc.get("okg_node_ids", {}).get("software_operational")),
            "okg_paper": (sc.get("okg_node_ids", {}).get("paper")
                           or sc.get("okg_node_ids", {}).get("paper_canonical")),
            "okg_ld_panel": sc.get("okg_node_ids", {}).get("ld_panel"),
        }
        # 1) Explicit key.
        for key in ("cs_path", "weights_path", "weight_path", "pip_file"):
            v = sc.get(key) or sc.get("paths", {}).get(key)
            if v and Path(v).exists():
                meta = {"source": "finemap_skill_sidecar",
                        "sidecar_path": str(sc_path),
                        "pip_file_key": key, **okg_passthrough}
                return Path(v), meta
        # 2) Auto-derive from output_prefix + sushie cs.tsv + gwas_extras.tsv.
        prefix = sc.get("output_prefix")
        if prefix:
            prefix_path = Path(prefix)
            cs_tsv = prefix_path.with_name(prefix_path.name + ".sushie.cs.tsv")
            extras_tsv = prefix_path.with_name(prefix_path.name + ".gwas_extras.tsv")
            if cs_tsv.exists() and extras_tsv.exists():
                joined = _materialize_rsid_pip(cs_tsv, extras_tsv, work_dir)
                meta = {"source": "finemap_skill_sidecar_autojoin",
                        "sidecar_path": str(sc_path),
                        "cs_tsv": str(cs_tsv),
                        "gwas_extras_tsv": str(extras_tsv),
                        **okg_passthrough}
                return joined, meta
            # 3) Fallback: weights.tsv with sushie_pip_all (already-rsid case)
            weights_tsv = prefix_path.with_name(prefix_path.name + ".sushie.weights.tsv")
            if weights_tsv.exists():
                meta = {"source": "finemap_skill_sidecar_weights",
                        "sidecar_path": str(sc_path),
                        "weights_tsv": str(weights_tsv),
                        **okg_passthrough}
                return weights_tsv, meta
        sys.exit(f"--finemap-sidecar {sc_path}: could not locate a PIP file. "
                 f"Tried explicit keys (cs_path/weights_path/weight_path/pip_file) "
                 f"and output_prefix-derived (<prefix>.sushie.cs.tsv + "
                 f"<prefix>.gwas_extras.tsv | <prefix>.sushie.weights.tsv).")
    return None, meta


def okg_lookup_dataset(okg_repo: Path, dataset_id: str) -> dict | None:
    """Best-effort OKG MCP query; returns attrs dict on hit, None otherwise."""
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
                         "clientInfo": {"name": "locuszoom", "version": "0.1"}}})
        read()
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
              "params": {"name": "get_node", "arguments": {"node_id": dataset_id}}})
        resp = read()
    finally:
        try: proc.stdin.close()
        except Exception: pass
        try: proc.wait(timeout=3)
        except subprocess.TimeoutExpired: proc.terminate()
    sc = resp.get("result", {}).get("structuredContent") or {}
    return sc.get("attrs") or None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sumstats", type=Path, required=True)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--gene", type=str)
    group.add_argument("--region", type=str)
    group.add_argument("--lead-snp", type=str)
    p.add_argument("--flank", type=int, default=100_000)
    p.add_argument("--ld-pop", type=str, default="EUR",
                   choices=["EUR", "AFR", "EAS", "AMR", "SAS"])
    p.add_argument("--build", type=str, default="hg19", choices=["hg19", "hg38"])
    p.add_argument("--ldlink-token", type=str,
                   help="LDlink API token (cached to ~/.cache/locuszoom/ldlink_token on first use)")
    p.add_argument("--out", type=Path, required=True,
                   help="Output prefix; writes <prefix>.pdf, .png, .locuszoom.json")
    p.add_argument("--okg-dataset-id", type=str,
                   help="OKG dataset_metadata node ID; recorded in sidecar.")
    p.add_argument("--okg-repo", type=Path,
                   default=Path(os.environ["OKG_REPO"])
                           if os.environ.get("OKG_REPO") else None)
    p.add_argument("--force-sparse", action="store_true",
                   help="Allow loci with fewer than 50 SNPs in window.")
    p.add_argument("--finemap-pip", type=Path,
                   help="TSV with per-variant PIP (and optional CS) columns. "
                        "Adds a PIP scatter panel; re-anchors lead to highest-PIP variant. "
                        "Column auto-detection: snp/rsid/SNP, pip/PIP/pip_all/sushie_pip_all, "
                        "cs/CSIndex/credible_set.")
    p.add_argument("--finemap-sidecar", type=Path,
                   help="A `.finemap.json` sidecar from the `finemap` skill; the PIP TSV "
                        "is derived from its `cs_path` or `weights_path` field.")
    args = p.parse_args()

    if not (args.gene or args.region or args.lead_snp):
        sys.stderr.write("BAD_LOCUS: pass one of --gene / --region / --lead-snp\n")
        sys.exit(EXIT_BAD_LOCUS)

    token, token_source = resolve_token(args.ldlink_token)
    sys.stderr.write(f"LDlink token source: {token_source}\n")

    args.sumstats = args.sumstats.expanduser().resolve()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if not args.sumstats.exists():
        sys.exit(f"sumstats not found: {args.sumstats}")

    work_dir = Path(tempfile.mkdtemp(prefix="locuszoom_"))
    norm_path = normalise_sumstats(args.sumstats, work_dir)

    pip_path, finemap_meta = resolve_finemap_pip(args, work_dir)
    if pip_path:
        sys.stderr.write(f"finemap PIP file: {pip_path}\n")

    region_args: list[str] = []
    locus_mode = ""
    if args.gene:
        locus_mode = "gene"
        region_args = ["--gene", args.gene, "--flank", str(args.flank)]
    elif args.region:
        locus_mode = "region"
        chrom, start, end = parse_region(args.region)
        region_args = ["--chrom", chrom, "--start", str(start), "--end", str(end)]
    elif args.lead_snp:
        locus_mode = "lead_snp"
        region_args = ["--lead-snp", args.lead_snp, "--flank", str(args.flank)]

    r_script = HERE / "locuszoom_run.R"
    cmd = ["Rscript", str(r_script),
           "--sumstats", str(norm_path),
           "--build", args.build,
           "--ld-pop", args.ld_pop,
           "--ldlink-token", token,
           "--out", str(args.out),
           "--force-sparse", "1" if args.force_sparse else "0",
           *region_args]
    if pip_path:
        cmd += ["--finemap-pip", str(pip_path)]
    sys.stderr.write(f"running R: {' '.join(cmd[:4])} ... (token redacted)\n")
    r_summary_path = work_dir / "r_summary.json"
    env = os.environ.copy()
    env["LOCUSZOOM_R_SUMMARY"] = str(r_summary_path)
    rc = subprocess.call(cmd, env=env)
    if rc == EXIT_TOKEN_INVALID:
        sys.stderr.write(
            "LDLINK_TOKEN_INVALID\n"
            "LDlink rejected the cached token (HTTP 401/403). Re-invoke with\n"
            "`--ldlink-token <fresh_token>` to refresh, or delete\n"
            f"{TOKEN_FILE} and retry.\n"
        )
        sys.exit(EXIT_TOKEN_INVALID)
    if rc != 0:
        sys.exit(f"R worker failed with exit {rc}; see stderr above")

    r_summary = (json.loads(r_summary_path.read_text())
                 if r_summary_path.exists() else {})

    okg_attrs = None
    if args.okg_repo and args.okg_dataset_id:
        okg_attrs = okg_lookup_dataset(args.okg_repo, args.okg_dataset_id)

    okg_node_ids = {"software": "software:locuszoomr"}
    if args.okg_dataset_id:
        okg_node_ids["dataset"] = args.okg_dataset_id
    okg_node_ids["ld_resource"] = f"external:ldlink_1000g_{args.ld_pop}"
    # Propagate any OKG nodes carried by the finemap sidecar.
    for k in ("okg_method", "okg_software", "okg_paper", "okg_ld_panel"):
        v = finemap_meta.get(k)
        if v:
            okg_node_ids[k.replace("okg_", "finemap_")] = v

    finemap_block = None
    fm_summary = r_summary.get("finemap") or {}
    if pip_path:
        finemap_block = {
            "pip_file": str(pip_path),
            "pip_file_sha256": sha256(pip_path),
            **finemap_meta,
            "n_variants_with_pip": fm_summary.get("n_variants_with_pip"),
            "max_pip": fm_summary.get("max_pip"),
            "pip_lead_snp": fm_summary.get("pip_lead_snp"),
            "n_credible_sets": fm_summary.get("n_credible_sets"),
            "credible_set_sizes": fm_summary.get("credible_set_sizes"),
        }

    manifest = {
        "sumstats_input": str(args.sumstats),
        "sumstats_sha256": sha256(args.sumstats),
        "output_pdf": f"{args.out}.pdf",
        "output_png": f"{args.out}.png",
        "locus_selection": {
            "mode": locus_mode,
            "gene": args.gene,
            "region": args.region,
            "lead_snp": args.lead_snp,
            "flank_bp": args.flank,
            **(r_summary.get("locus") or {}),
        },
        "build": args.build,
        "ensembl_db": ("EnsDb.Hsapiens.v75" if args.build == "hg19"
                        else "EnsDb.Hsapiens.v86"),
        "ld": {
            "source": "LDlink",
            "endpoint": "https://ldlink.nih.gov/LDlinkRest/ldproxy",
            "population": args.ld_pop,
            "lead_snp": r_summary.get("lead_snp"),
            "n_ld_pairs": r_summary.get("n_ld_pairs"),
        },
        "n_snps_in_window": r_summary.get("n_snps_in_window"),
        "finemap": finemap_block,
        "okg_node_ids": okg_node_ids,
        "okg_dataset_attrs": okg_attrs,
        "locuszoomr_version": r_summary.get("locuszoomr_version"),
        "captured_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    sidecar = Path(f"{args.out}.locuszoom.json")
    sidecar.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"plot -> {args.out}.pdf and {args.out}.png")
    print(f"sidecar -> {sidecar}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
