#!/usr/bin/env python3
"""
snp_report.py — extract common + thyroid-related SNPs from consumer raw DNA data
(23andMe or AncestryDNA .txt / .zip) and build an HTML report plus CSV extracts.

Usage
-----
    python snp_report.py                       # auto-finds genome_*.txt / .zip next to this folder
    python snp_report.py path/to/genome.txt    # or .zip
    python snp_report.py genome.txt -o out_dir --offline

Outputs (in the output folder, default = this script's folder)
    snp_report.html                 human-readable report
    snp_extracted.csv               every curated SNP with your genotype + interpretation
    thyroid_gene_region_snps.csv    every SNP on your chip inside thyroid-pathway genes
                                    (annotated or not)

The SNP list and interpretations live in snp_database.py — edit that file to add SNPs.
Gene coordinates for the region export come from Ensembl (GRCh37) and are cached in
gene_coords_grch37.json; only gene symbols are sent, never genotype data.
"""

import argparse
import csv
import gzip
import datetime as dt
import html
import io
import json
import re
import sys
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from snp_database import (APOE_NOTES, APOE_TABLE, CATEGORIES, SNPS,
                          THYROID_GENES)

HERE = Path(__file__).resolve().parent
GENE_CACHE = HERE / "gene_coords_grch37.json"
GENE_FLANK = 5_000  # bp added each side of a gene for the region export
POS_CACHE = HERE / "snp_positions_grch37.json"
NO_CALLS = {"--", "00", "", "-", "0"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CHR = {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT", "XY": "X", "PAR": "X"}
BUILD_ANCHOR = {742429: "NCBI36", 752566: "GRCh37", 817186: "GRCh38"}  # position of rs3094315
HEAD = {
    "id": re.compile(r"^(rsid|rs id|snp name|snp|name|marker|id)$"),
    "chr": re.compile(r"^(chr|chromosome|chrom)$"),
    "pos": re.compile(r"^(position|pos|coordinate)$"),
    "gt": re.compile(r"^(genotype|result|call)$"),
    "a1": re.compile(r"^allele ?_?1"),
    "a2": re.compile(r"^allele ?_?2"),
}


def is_indel(gt):
    return bool(gt) and set(gt) <= {"I", "D"}


# --------------------------------------------------------------------------- input
def find_genome_file():
    for base in (Path.cwd(), HERE.parent, HERE):
        for pattern in ("**/genome*.txt", "**/*23andme*.txt", "**/AncestryDNA*.txt",
                        "genome*.zip", "*AncestryDNA*.zip"):
            hits = sorted(base.glob(pattern))
            if hits:
                return hits[0]
    return None


def open_raw(path: Path):
    """Yield text lines from a .txt/.csv, a .gz (all members), or the data file inside a .zip."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"PK":
        with zipfile.ZipFile(path) as zf:
            names = [n for n in zf.namelist()
                     if re.search(r"\.(txt|csv|tsv)$", n, re.I) and "__MACOSX" not in n]
            if not names:
                sys.exit("No .txt or .csv file found inside the zip.")
            gfg = [n for n in names if re.search(r"GFG.*23andMe\.txt$", n, re.I)]
            name = gfg[0] if gfg else max(names, key=lambda n: zf.getinfo(n).file_size)
            with zf.open(name) as fh:
                yield from io.TextIOWrapper(fh, encoding="utf-8-sig", errors="replace")
    elif magic == b"\x1f\x8b":
        with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as fh:  # reads every member
            yield from fh
    else:
        with open(path, encoding="utf-8-sig", errors="replace") as fh:
            yield from fh


def split_fields(line):
    line = line.replace('"', "")
    parts = re.split(r"\s*[\t,]\s*", line) if re.search(r"[\t,]", line) else line.split()
    return [x.strip() for x in parts if x.strip()]


class LineParser:
    """Turns raw-data lines into (rsid, chrom, pos, genotype). Handles 23andMe, AncestryDNA,
    MyHeritage, FamilyTreeDNA, Living DNA, Sano and other header-row formats."""

    def __init__(self):
        self.cols = None
        self.build = "unknown"

    def __call__(self, line):
        line = line.strip()
        if not line:
            return None
        if line.startswith("##fileformat=VCF") or line.startswith("#CHROM\tPOS\tID"):
            sys.exit("This is a VCF sequencing file; only genotyping-chip raw data files are supported.")
        if line[0] == "#":
            low = line.lower()
            if re.search(r"build 36|ncbi36|hg18", low):
                self.build = "NCBI36"
            elif re.search(r"build 37|grch37|hg19", low):
                self.build = "GRCh37"
            elif re.search(r"build 38|grch38|hg38", low):
                self.build = "GRCh38"
            return None
        f = split_fields(line)
        if len(f) < 3:
            return None
        pos_i = self.cols["pos"] if self.cols else 2
        if pos_i >= len(f) or not f[pos_i].isdigit():          # header row: learn its columns
            low = [x.lower() for x in f]
            idx = {k: next((i for i, x in enumerate(low) if rx.match(x)), -1) for k, rx in HEAD.items()}
            if min(idx["id"], idx["chr"], idx["pos"]) >= 0 and max(idx["gt"], idx["a1"]) >= 0:
                self.cols = idx
            return None
        try:
            if self.cols:
                c = self.cols
                rsid, chrom, pos = f[c["id"]], f[c["chr"]], int(f[c["pos"]])
                gt = f[c["gt"]] if c["gt"] >= 0 else f[c["a1"]] + (f[c["a2"]] if c["a2"] < len(f) else "")
            else:
                rsid, chrom, pos = f[0], f[1], int(f[2])
                if len(f) >= 5 and len(f[3]) <= 1 and len(f[4]) <= 1:
                    gt = f[3] + f[4]
                else:
                    gt = f[3] if len(f) > 3 else ""
        except (IndexError, ValueError):
            return None
        gt = re.sub(r"[0\-]", "", gt.upper()) or "--"
        chrom = chrom.upper()
        chrom = chrom[3:] if chrom.startswith("CHR") else chrom
        chrom = CHR.get(chrom, chrom)
        if chrom == "0":
            return None
        return rsid, chrom, pos, gt


# --------------------------------------------------------------------------- gene regions
def load_gene_coords(genes, offline=False):
    cache = {}
    if GENE_CACHE.exists():
        cache = json.loads(GENE_CACHE.read_text(encoding="utf-8"))
    missing = [g for g in genes if g not in cache]
    if missing and not offline:
        try:
            req = urllib.request.Request(
                "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens",
                data=json.dumps({"symbols": missing}).encode(),
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
            for sym, rec in data.items():
                if rec and rec.get("seq_region_name"):
                    cache[sym] = {"chrom": rec["seq_region_name"], "start": rec["start"],
                                  "end": rec["end"], "strand": rec.get("strand")}
            GENE_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True), encoding="utf-8")
        except Exception as exc:  # network failure is non-fatal
            print(f"  ! could not fetch gene coordinates from Ensembl ({exc}); "
                  f"region export limited to cached genes", file=sys.stderr)
    return {g: cache[g] for g in genes if g in cache}


# --------------------------------------------------------------------------- scoring
def interpret(snp, gt):
    """Return dict with copies, status, label, note for one SNP."""
    if gt is None:
        return dict(copies="", status="missing", label="Not on chip", note="")
    if gt in NO_CALLS:
        return dict(copies="", status="nocall", label="No call", note="")
    if is_indel(gt):
        return dict(copies="", status="nocall", label="Indel call",
                    note="Your file reports an insertion/deletion call here, not a SNP genotype.")
    expected = set(snp["alleles"])
    raw, flipped = gt, False
    if not set(gt) <= expected:
        comp = gt.translate(COMPLEMENT)
        palindromic = snp["alleles"][0].translate(COMPLEMENT) == snp["alleles"][1]
        if palindromic or not set(comp) <= expected:
            return dict(copies="", status="check", label="Check strand",
                        note=f"Genotype {gt} doesn't match expected alleles "
                             f"{'/'.join(snp['alleles'])}; possible strand flip or different variant.")
        gt, flipped = comp, True
    ploidy = len(gt)
    n = gt.count(snp["effect"])
    frac = n / ploidy
    note = snp["notes"].get(n if ploidy == 2 else n * 2, "")
    if flipped:
        note = (note + " " if note else "") +             f"(Your file reported {raw} on the opposite strand; shown here as {gt}.)"
    kind = snp["kind"]
    if ploidy == 1:
        copy_txt = f"{n} (hemizygous)"
    else:
        copy_txt = str(n)
    if kind == "risk":
        status = "typical" if n == 0 else ("hom" if frac == 1 else "het")
        label = {"typical": "Typical", "het": "1 copy", "hom": "2 copies"}[status]
        if ploidy == 1 and n:
            label = "1 copy (X-linked)"
    elif kind == "benefit":
        status = "good" if n else "neutral"
        label = "Favourable" if n else "Typical"
    else:
        status = "info"
        label = "Info"
    return dict(gt=gt, copies=copy_txt, status=status, label=label, note=note)


def apoe_type(found):
    a = found.get("rs429358")
    b = found.get("rs7412")
    if not a or not b or a in NO_CALLS or b in NO_CALLS:
        return None
    key = ("".join(sorted(a)), "".join(sorted(b)))
    return APOE_TABLE.get(key, f"Unrecognised combination ({a}, {b})")


# --------------------------------------------------------------------------- main pass
def scan(path, regions):
    wanted = {s["rsid"] for s in SNPS}
    positions = json.loads(POS_CACHE.read_text(encoding="utf-8")) if POS_CACHE.exists() else {}
    by_pos = {(c, p): sid for sid, (c, p) in positions.items()}
    found, pos_of, via_pos = {}, {}, set()
    region_hits = []
    by_chrom = defaultdict(list)
    for gene, r in regions.items():
        by_chrom[r["chrom"]].append((r["start"] - GENE_FLANK, r["end"] + GENE_FLANK, gene))
    stats = Counter()
    parse = LineParser()
    anchor = None
    for line in open_raw(path):
        rec = parse(line)
        if not rec:
            continue
        rsid, chrom, pos, gt = rec
        stats["total"] += 1
        if gt in NO_CALLS:
            stats["nocall"] += 1
        if rsid == "rs3094315":
            anchor = pos
        if rsid in wanted:
            if rsid not in found or rsid in via_pos:
                found[rsid], pos_of[rsid] = gt, (chrom, pos)
                via_pos.discard(rsid)
        else:  # provider-specific ID (or other rsID) for a report SNP: match by GRCh37 position
            sid = by_pos.get((chrom, pos))
            if sid and sid not in found:
                found[sid], pos_of[sid] = gt, (chrom, pos)
                via_pos.add(sid)
        for start, end, gene in by_chrom.get(chrom, ()):
            if start <= pos <= end:
                region_hits.append((gene, rsid, chrom, pos, gt))
    build = parse.build
    if build == "unknown" and anchor:
        build = BUILD_ANCHOR.get(anchor, "unknown")
    if build not in ("GRCh37", "unknown"):  # positions only valid on GRCh37
        for sid in via_pos:
            found.pop(sid, None)
            pos_of.pop(sid, None)
        region_hits = []
    return found, pos_of, region_hits, stats, build


# --------------------------------------------------------------------------- outputs
def write_extract_csv(rows, out):
    cols = ["group", "section", "gene", "rsid", "variant", "chrom", "position", "genotype",
            "expected_alleles", "effect_allele", "effect_copies", "status", "evidence",
            "interpretation", "genotype_note"]
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([r["group"], r["section"], r["gene"], r["rsid"], r["variant"], r["chrom"],
                        r["pos"], r["gt"] or "", "/".join(r["alleles"]), r["effect"], r["copies"],
                        r["label"], r["evidence"], r["desc"], r["note"]])


def write_region_csv(region_hits, curated, out):
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "rsid", "chrom", "position", "genotype", "in_curated_report"])
        for gene, rsid, chrom, pos, gt in sorted(region_hits, key=lambda x: (x[0], x[3])):
            w.writerow([gene, rsid, chrom, pos, gt, "yes" if rsid in curated else ""])


CSS = """
:root{--bg:#f7f7f5;--card:#fff;--ink:#1d1d1f;--muted:#6b6b70;--line:#e4e4e0;--accent:#2f5d8a;
--typ-bg:#e9f5ec;--typ:#1f6b37;--het-bg:#fff4dc;--het:#8a5a00;--hom-bg:#fde8e6;--hom:#a3261b;
--good-bg:#e4f1fb;--good:#1d5d8f;--info-bg:#efeff3;--info:#55555c;--chk-bg:#f3e8fb;--chk:#6c2f91}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#141416;--card:#1d1d20;
--ink:#ececef;--muted:#9a9aa2;--line:#2e2e33;--accent:#8db8e0;--typ-bg:#183322;--typ:#86d29d;
--het-bg:#3a2d10;--het:#f2c46b;--hom-bg:#3d1a17;--hom:#ff9b90;--good-bg:#16283a;--good:#8fc4ef;
--info-bg:#26262b;--info:#b4b4bb;--chk-bg:#2e1f3a;--chk:#d3a6f2}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 -apple-system,"Segoe UI",
Roboto,Helvetica,Arial,sans-serif}
main{max-width:1180px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:28px;margin:0 0 4px}h2{font-size:21px;margin:40px 0 12px;padding-top:8px;
border-top:2px solid var(--line)}h3{font-size:16px;margin:28px 0 8px;color:var(--accent)}
.sub{color:var(--muted);margin:0 0 20px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:12px;margin:16px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card>b{display:block;font-size:13px;color:var(--muted);font-weight:600;text-transform:uppercase;
letter-spacing:.03em}.card .big{font-size:22px;font-weight:700;margin-top:2px}
.warn{background:var(--het-bg);border:1px solid var(--line);border-radius:10px;padding:12px 16px;
font-size:14px}
ul.find{margin:6px 0 0;padding-left:18px}ul.find li{margin:3px 0}
.tbl{overflow-x:auto;background:var(--card);border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:top;text-align:left}
th{font-size:12px;text-transform:uppercase;letter-spacing:.03em;color:var(--muted);
background:var(--card);position:sticky;top:0}
tr:last-child td{border-bottom:0}
td.gt{font-family:Consolas,Menlo,monospace;font-weight:700;font-size:15px;white-space:nowrap}
td.rs{font-family:Consolas,Menlo,monospace;font-size:12.5px;white-space:nowrap}
td.rs a{color:var(--accent);text-decoration:none}
.chip{display:inline-block;padding:2px 8px;border-radius:99px;font-size:12px;font-weight:600;
white-space:nowrap}
.typical,.neutral{background:var(--typ-bg);color:var(--typ)}.het{background:var(--het-bg);color:var(--het)}
.hom{background:var(--hom-bg);color:var(--hom)}.good{background:var(--good-bg);color:var(--good)}
.info,.missing,.nocall{background:var(--info-bg);color:var(--info)}.check{background:var(--chk-bg);color:var(--chk)}
.ev{font-size:12px;color:var(--muted)}.note{display:block;margin-top:4px;font-weight:600}
.desc{min-width:320px}
details{margin-top:8px}summary{cursor:pointer;color:var(--accent)}
.small{font-size:13px;color:var(--muted)}
nav{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:14px;margin:8px 0 0}
nav a{color:var(--accent)}
.sub,p,li{overflow-wrap:anywhere}
h2,h3{scroll-margin-top:56px}
td .lbl{display:none}
@media (max-width:720px){
 body{font-size:15px}
 main{padding:16px 16px 48px}
 h1{font-size:22px;line-height:1.25}
 h2{font-size:19px;margin-top:32px}
 nav{position:sticky;top:0;z-index:5;flex-wrap:nowrap;overflow-x:auto;gap:8px;margin:12px -16px 0;
  padding:8px 16px;background:var(--bg);border-bottom:1px solid var(--line);
  scrollbar-width:none;-webkit-overflow-scrolling:touch}
 nav::-webkit-scrollbar{display:none}
 nav a{flex:0 0 auto;white-space:nowrap;text-decoration:none;font-size:13px;padding:6px 12px;
  border-radius:99px;background:var(--card);border:1px solid var(--line)}
 .cards{grid-template-columns:1fr 1fr;gap:10px}
 .cards .card.wide{grid-column:1/-1}
 .card{padding:12px 14px}.card .big{font-size:20px}
 ul.find{padding-left:0;list-style:none}
 ul.find li{margin:0;padding:10px 0;border-top:1px solid var(--line)}
 ul.find li:first-child{border-top:0}
 ul.find .chip{margin-bottom:4px}
 /* SNP tables become one card per row */
 .tbl.snps{background:none;border:0;overflow:visible}
 .tbl.snps table,.tbl.snps tbody{display:block}
 .tbl.snps thead{display:none}
 .tbl.snps tr{display:grid;grid-template-columns:1fr 1fr auto;column-gap:12px;row-gap:6px;
  grid-template-areas:"g g res" "v v res" "rs you eff" "desc desc desc";
  background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin:0 0 10px}
 .tbl.snps td{display:block;padding:0;border:0;min-width:0}
 .tbl.snps td.g{grid-area:g;font-size:16px}
 .tbl.snps td.v{grid-area:v;color:var(--muted);font-size:13px;margin-top:-4px}
 .tbl.snps td.res{grid-area:res;text-align:right}
 .tbl.snps td.rs{grid-area:rs}.tbl.snps td.you{grid-area:you}.tbl.snps td.eff{grid-area:eff;text-align:right}
 .tbl.snps td.rs,.tbl.snps td.you,.tbl.snps td.eff{background:var(--bg);border-radius:8px;padding:6px 8px}
 .tbl.snps td .lbl{display:block;font:600 10.5px/1.3 -apple-system,"Segoe UI",Roboto,sans-serif;
  text-transform:uppercase;letter-spacing:.04em;color:var(--muted);margin-bottom:2px}
 .tbl.snps td.eff .ev{display:none}
 .tbl.snps td.desc{grid-area:desc;min-width:0;font-size:14px;line-height:1.5}
 .tbl.snps td.rs a{display:inline-block;padding:2px 0}
 td.rs{white-space:normal}
}
@media (max-width:380px){.cards{grid-template-columns:1fr}
 .tbl.snps tr{grid-template-areas:"g g res" "v v v" "rs you eff" "desc desc desc"}}
@media print{body{background:#fff}.tbl{overflow:visible}th{position:static}h2{break-before:page}}
"""


def first_sentence(text):
    m = re.search(r"(?<!i\.e)(?<!e\.g)\.\s+(?=[A-Z])", text)
    return text[:m.start() + 1] if m else text


def chip(status, label):
    return f'<span class="chip {status}">{html.escape(label)}</span>'


def build_html(rows, apoe, region_hits, region_genes, stats, build, src, out):
    e = html.escape
    sec_rows = defaultdict(list)
    for r in rows:
        sec_rows[r["cat"]].append(r)

    present = [r for r in rows if r["status"] not in ("missing",)]
    flagged = [r for r in rows if r["status"] in ("het", "hom") and r["evidence"] != "limited"]
    favourable = [r for r in rows if r["status"] == "good"]
    checks = [r for r in rows if r["status"] == "check"]

    def finding_li(r):
        txt = f'<b>{e(r["gene"])} {e(r["variant"])}</b> ({e(r["rsid"])}, {e(r["gt"])}) — '
        txt += e(r["note"] or first_sentence(r["desc"]))
        return f"<li>{chip(r['status'], r['label'])} {txt}</li>"

    parts = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
             f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>SNP Report</title><style>{CSS}</style></head><body><main>"]
    parts.append("<h1>Genetic SNP Report — Thyroid &amp; Common Variants</h1>")
    parts.append(f"<p class='sub'>Source: {e(src.name)} · generated {dt.date.today():%d %b %Y} · "
                 f"reference build {e(build)}</p>")

    nav = " ".join(f"<a href='#{k}'>{e(t)}</a>" for k, t, _ in CATEGORIES)
    parts.append(f"<nav><a href='#summary'>Summary</a> {nav} <a href='#regions'>All thyroid-gene SNPs</a>"
                 f" <a href='#about'>Method &amp; caveats</a></nav>")

    # ---- summary
    parts.append("<h2 id='summary'>Summary</h2><div class='cards'>")
    parts.append(f"<div class='card'><b>SNPs in raw file</b><div class='big'>{stats['total']:,}</div>"
                 f"<span class='small'>no-call rate {stats['nocall'] / max(stats['total'], 1):.2%}</span></div>")
    parts.append(f"<div class='card'><b>Curated SNPs found</b><div class='big'>{len(present)} / {len(rows)}</div>"
                 f"<span class='small'>others aren't on this chip version</span></div>")
    parts.append(f"<div class='card'><b>Thyroid-gene SNPs on chip</b><div class='big'>{len(region_hits):,}</div>"
                 f"<span class='small'>across {len(region_genes)} genes (see CSV)</span></div>")
    if apoe:
        parts.append(f"<div class='card wide'><b>APOE type</b><div class='big'>{e(apoe)}</div>"
                     f"<span class='small'>{e(APOE_NOTES.get(apoe, ''))}</span></div>")
    parts.append("</div>")

    # thyroid at a glance
    parts.append("<h3>Thyroid at a glance</h3><div class='card'>")
    for key, title, group in CATEGORIES:
        if group != "Thyroid":
            continue
        rs_ = sec_rows[key]
        hits = [r for r in rs_ if r["status"] in ("het", "hom", "good")]
        parts.append(f"<p style='margin:10px 0 2px'><b>{e(title)}</b> — "
                     f"{sum(r['status'] != 'missing' for r in rs_)} SNPs tested</p>")
        if hits:
            parts.append("<ul class='find'>" + "".join(finding_li(r) for r in hits) + "</ul>")
        else:
            parts.append("<p class='small' style='margin:0'>No variant alleles flagged.</p>")
    parts.append("</div>")

    common_flags = [r for r in flagged if r["group"] == "Common"]
    if common_flags or favourable:
        parts.append("<h3>Other notable findings (moderate/strong evidence)</h3><div class='card'><ul class='find'>")
        parts.extend(finding_li(r) for r in common_flags + [f for f in favourable if f["group"] == "Common"])
        parts.append("</ul></div>")
    if checks:
        parts.append("<h3>Needs checking</h3><div class='card'><ul class='find'>")
        parts.extend(f"<li>{chip('check', 'Check')} {e(r['gene'])} {e(r['rsid'])}: {e(r['note'])}</li>"
                     for r in checks)
        parts.append("</ul></div>")

    parts.append("<p class='warn' style='margin-top:16px'><b>How to read this:</b> "
                 "<span class='chip typical'>Typical</span> no copies of the effect allele · "
                 "<span class='chip het'>1 copy</span> heterozygous · "
                 "<span class='chip hom'>2 copies</span> homozygous · "
                 "<span class='chip good'>Favourable</span> carries a protective/beneficial allele · "
                 "<span class='chip info'>Info</span> direction unclear or descriptive. "
                 "Most common variants shift risk by small amounts; they are not diagnoses. "
                 "Your lab results (TSH, FT4, FT3, rT3, antibodies) always outweigh these SNPs.</p>")

    # ---- sections
    current_group = None
    for key, title, group in CATEGORIES:
        if group != current_group:
            parts.append(f"<h2>{e(group)} variants</h2>")
            current_group = group
        parts.append(f"<h3 id='{key}'>{e(title)}</h3><div class='tbl snps'><table><thead><tr>"
                     "<th>Gene</th><th>SNP</th><th>Variant</th><th>Your genotype</th>"
                     "<th>Effect allele</th><th>Result</th><th>What it means</th></tr></thead><tbody>")
        for r in sec_rows[key]:
            link = f"https://www.snpedia.com/index.php/{r['rsid']}" if r["rsid"].startswith("rs") else ""
            rs_html = f"<a href='{link}' target='_blank' rel='noopener'>{e(r['rsid'])}</a>" if link else e(r["rsid"])
            note = f"<span class='note'>{e(r['note'])}</span>" if r["note"] else ""
            parts.append(
                f"<tr><td class='g'><b>{e(r['gene'])}</b></td>"
                f"<td class='rs'><span class='lbl'>SNP</span>{rs_html}</td><td class='v'>{'' if r['variant'] == r['rsid'] else e(r['variant'])}</td>"
                f"<td class='gt you'><span class='lbl'>You</span>{e(r['gt'] or '—')}</td>"
                f"<td class='gt eff'><span class='lbl'>Effect</span>{e(r['effect'])}"
                f"<div class='ev'>{e('/'.join(r['alleles']))}</div></td>"
                f"<td class='res'>{chip(r['status'], r['label'])}<div class='ev'>{e(r['evidence'])} evidence</div></td>"
                f"<td class='desc'>{e(r['desc'])}{note}</td></tr>")
        parts.append("</tbody></table></div>")

    # ---- regions
    parts.append("<h2 id='regions'>All thyroid-pathway gene SNPs on your chip</h2>")
    parts.append("<p>Every SNP your chip measured inside these genes (± 5 kb) is exported to "
                 "<code>thyroid_gene_region_snps.csv</code>. Most have no well-established effect — "
                 "the table above contains the ones with published associations.</p>")
    per_gene = Counter(h[0] for h in region_hits)
    parts.append("<div class='tbl'><table><thead><tr><th>Gene</th><th>Location (GRCh37)</th>"
                 "<th>SNPs on chip</th></tr></thead><tbody>")
    for g in sorted(region_genes):
        r = region_genes[g]
        parts.append(f"<tr><td><b>{e(g)}</b></td><td class='rs'>chr{e(str(r['chrom']))}:"
                     f"{r['start']:,}-{r['end']:,}</td><td>{per_gene.get(g, 0)}</td></tr>")
    parts.append("</tbody></table></div>")

    # ---- about
    parts.append("<h2 id='about'>Method &amp; caveats</h2><div class='card small'>"
                 "<p>Genotypes are read from the raw file on the GRCh37 plus strand and compared with a "
                 "curated list (snp_database.py). For each SNP the report counts copies of the "
                 "<i>effect allele</i> named in the table. If a genotype contains alleles other than those "
                 "expected, it is flagged rather than interpreted.</p>"
                 "<p>Consumer arrays are research-grade: rare variants can be mis-called, copy-number "
                 "changes and deletions (e.g. GSTM1/GSTT1, CYP2D6 duplications) are not detected, and "
                 "many SNPs are missing from older chip versions. Confirm anything clinically important "
                 "(e.g. APOE, Factor V Leiden, HFE, pharmacogenomic results) with a clinical-grade test "
                 "before acting on it.</p>"
                 "<p>Effect sizes for thyroid SNPs are small: DIO1/DIO2 variants shift average hormone "
                 "levels by a few percent, and autoimmunity SNPs raise relative, not absolute, risk. "
                 "This report is informational and is not medical advice.</p></div>")
    parts.append("</main></body></html>")
    out.write_text("\n".join(parts), encoding="utf-8")


# --------------------------------------------------------------------------- cli
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome", nargs="?", help="raw data file (.txt or .zip)")
    ap.add_argument("-o", "--outdir", default=str(HERE), help="output folder (default: script folder)")
    ap.add_argument("--offline", action="store_true", help="don't contact Ensembl for gene coordinates")
    ap.add_argument("--no-gene-regions", action="store_true", help="skip the thyroid gene-region export")
    args = ap.parse_args()
    for stream in (sys.stdout, sys.stderr):  # Windows consoles default to cp1252
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")

    src = Path(args.genome) if args.genome else find_genome_file()
    if not src or not src.exists():
        sys.exit("No genome file found. Pass the path to your raw data .txt or .zip.")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {src} ...")
    regions = {} if args.no_gene_regions else load_gene_coords(THYROID_GENES, args.offline)
    found, pos_of, region_hits, stats, build = scan(src, regions)
    if build not in ("GRCh37", "unknown"):
        print(f"  ! file is {build}: rsID lookups still work, gene-region export skipped", file=sys.stderr)

    cat_meta = {k: (t, g) for k, t, g in CATEGORIES}
    rows = []
    for snp in SNPS:
        gt = found.get(snp["rsid"])
        res = interpret(snp, gt)
        chrom, pos = pos_of.get(snp["rsid"], ("", ""))
        title, group = cat_meta[snp["cat"]]
        rows.append({**snp, **res, "gt": res.get("gt", gt), "chrom": chrom, "pos": pos,
                     "section": title, "group": group})

    apoe = apoe_type({r["rsid"]: r["gt"] for r in rows})

    write_extract_csv(rows, outdir / "snp_extracted.csv")
    if regions:
        write_region_csv(region_hits, found, outdir / "thyroid_gene_region_snps.csv")
    build_html(rows, apoe, region_hits, regions, stats, build, src, outdir / "snp_report.html")

    n_found = sum(r["status"] != "missing" for r in rows)
    print(f"  {stats['total']:,} SNPs read; {n_found}/{len(rows)} curated SNPs present; "
          f"{len(region_hits):,} thyroid-gene SNPs")
    if apoe:
        print(f"  APOE: {apoe}")
    for r in rows:
        if r["status"] == "check":
            print(f"  ! {r['rsid']} {r['gene']}: {r['note']}")
    print(f"Wrote {outdir / 'snp_report.html'}")


if __name__ == "__main__":
    main()
