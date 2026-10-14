#!/usr/bin/env python3
"""
Build the static one-page site (docs/index.html) from template.html + the shared
SNP database, so the website and snp_report.py always use the same SNP list.

    python site/build_site.py                      # writes docs/index.html
    python site/build_site.py --preview out.html   # also writes a claude.ai Artifact preview

Settings come from site/config.json (kept out of git; copy config.example.json to
start). Everything you put there ends up visible in the public page source, so
it only holds public values:

  donate     Venmo / Cash App / PayPal. "@name", "$tag", "name" or a full https:// URL.
             Empty values are hidden.
  support    The personal card at the bottom: "title", "message" (list of paragraphs),
             "substack" URL, "linkedin" URL (shown as a small icon), "contact_url" + "contact_label" (e.g. a Tally or
             Google Form, so visitors can reach you without seeing your email), and "photo" (an image file in site/, embedded into the
             page at build time so the file itself stays out of git).
             With no donation links and no Substack, the card is omitted.
  analytics  {"provider": "cloudflare" | "goatcounter" | "", "id": "..."}
             cloudflare  -> id is the Web Analytics beacon token
             goatcounter -> id is your site code (the CODE in CODE.goatcounter.com)
             The page's Content-Security-Policy is opened for that one service only,
             and the privacy wording on the page is updated to mention it.
  repo_url   Link used for the "free, open-source" credit.
  site_url   The live site's address, printed at the top of PDF exports and used for the
             canonical link, social-preview tags and docs/sitemap.xml.
  verification  {"google": "...", "bing": "..."}: the content= values of the
             google-site-verification / msvalidate.01 meta tags that Google Search
             Console and Bing Webmaster Tools ask you to add. Empty values are skipped.

docs/ is what GitHub Pages serves (Settings > Pages > main /docs); the same
index.html also works on any other static host. There is no server component: all
analysis runs in the visitor's browser.
"""

import argparse
import base64
import html
import json
import re
import sys
import urllib.request
from pathlib import Path

SITE = Path(__file__).resolve().parent
ROOT = SITE.parent
sys.path.insert(0, str(ROOT))

from snp_database import APOE_NOTES, APOE_TABLE, CATEGORIES, SNPS, THYROID_GENES  # noqa: E402

DOCS = ROOT / "docs"
OG_IMAGE = "og-image.png"  # docs/og-image.png, the 1200x630 social-preview card
POS_CACHE = ROOT / "snp_positions_grch37.json"
# 23andMe internal IDs -> the public rsID used to look up their position
ID_ALIASES = {"i3002432": "rs1799963"}
DEFAULT_REPO = "https://github.com/aimoony/thyroid-snp-check"
DEFAULT_SITE = "https://aimoony.github.io/thyroid-snp-check/"

# Button icons (site/icons/*.svg) and the colour each brand mark is drawn in.
ICON_COLORS = {"cashapp": "#00C244", "paypal": "#009CDE", "substack": "#FF6719"}


def icon(name):
    path = SITE / "icons" / f"{name}.svg"
    if not path.exists():
        return ""
    svg = re.sub(r"<!--.*?-->|<title>.*?</title>", "", path.read_text(encoding="utf-8"), flags=re.S).strip()
    svg = re.sub(r"\s+", " ", svg)
    attrs = 'class="ico" aria-hidden="true" focusable="false"'
    if name in ICON_COLORS:
        attrs += f' fill="{ICON_COLORS[name]}"'
    head, rest = svg.split(">", 1)  # only touch the opening <svg> tag
    head = re.sub(r'\s(width|height|class|role)="[^"]*"', "", head)
    return head.replace("<svg", f"<svg {attrs}", 1) + ">" + rest


DONATE = {
    "venmo": ("Venmo", lambda v: f"https://venmo.com/u/{v.lstrip('@')}"),
    "cashapp": ("Cash App", lambda v: f"https://cash.app/${v.lstrip('$')}"),
    "paypal": ("PayPal", lambda v: f"https://www.paypal.me/{v}"),
}

ANALYTICS = {
    "cloudflare": dict(
        name="Cloudflare Web Analytics",
        tag='<script defer src="https://static.cloudflareinsights.com/beacon.min.js" '
            'data-cf-beacon=\'{{"token": "{id}"}}\'></script>',
        script="https://static.cloudflareinsights.com", connect="https://cloudflareinsights.com", img=""),
    "goatcounter": dict(
        name="GoatCounter",
        tag='<script data-goatcounter="https://{id}.goatcounter.com/count" async '
            'src="https://gc.zgo.at/count.js"></script>',
        script="https://gc.zgo.at", connect="https://{id}.goatcounter.com", img="https://{id}.goatcounter.com"),
}

# Privacy copy that must change when analytics is on, so the page never overclaims.
PRIVACY_WITH_ANALYTICS = [
    ("Analysed in your browser and never uploaded. The page is blocked from sending data anywhere.",
     "Analysed in your browser and never uploaded. The page can only send an anonymous page-view count to {name}, never your data."),
]


def load_config():
    for name in ("config.json", "config.example.json"):
        path = SITE / name
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8")), name
    return {}, None


def load_positions():
    """GRCh37 positions of report SNPs, so files that use their own SNP IDs still match."""
    cache = json.loads(POS_CACHE.read_text(encoding="utf-8")) if POS_CACHE.exists() else {}
    want = {s["rsid"]: ID_ALIASES.get(s["rsid"], s["rsid"]) for s in SNPS}
    missing = [rs for sid, rs in want.items() if sid not in cache]
    if missing:
        try:
            req = urllib.request.Request(
                "https://grch37.rest.ensembl.org/variation/homo_sapiens",
                data=json.dumps({"ids": missing}).encode(),
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read())
            back = {rs: sid for sid, rs in want.items()}
            for rs, rec in data.items():
                maps = [m for m in rec.get("mappings", [])
                        if m.get("assembly_name") == "GRCh37" and m.get("seq_region_name", "").isalnum()
                        and "_" not in m.get("seq_region_name", "")]
                if len(maps) == 1 and maps[0]["start"] == maps[0]["end"] and rs in back:
                    cache[back[rs]] = [maps[0]["seq_region_name"], maps[0]["start"]]
            POS_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True), encoding="utf-8")
        except Exception as exc:
            print(f"  ! couldn't fetch SNP positions from Ensembl ({exc}); position matching limited",
                  file=sys.stderr)
    return {sid: cache[sid] for sid in want if sid in cache}


def support_html(cfg):
    """The personal "support me" card at the bottom of the page: optional photo,
    message, donation buttons and a Substack link."""
    sup = cfg.get("support", {})
    links = []
    for key, (label, to_url) in DONATE.items():
        val = (cfg.get("donate", {}).get(key) or "").strip()
        if not val:
            continue
        url = val if val.startswith("https://") else to_url(val)
        links.append(f'<a class="btn ghost" href="{html.escape(url)}" target="_blank" '
                     f'rel="noopener noreferrer">{icon(key)}{label}</a>')
    if sup.get("contact_url"):
        links.insert(0, f'<a class="btn" href="{html.escape(sup["contact_url"])}" target="_blank" '
                     f'rel="noopener noreferrer">{icon("work")}{html.escape(sup.get("contact_label") or "Work with me")}</a>')
    if sup.get("substack"):
        cls = "btn ghost" if sup.get("contact_url") else "btn"
        links.append(f'<a class="{cls}" href="{html.escape(sup["substack"])}" target="_blank" '
                     f'rel="noopener noreferrer">{icon("substack")}Substack</a>')
    if not links:
        return ""
    photo = ""
    if sup.get("photo"):
        # Embedded as a data URI so the photo file itself never has to be committed.
        path = SITE / sup["photo"]
        if path.exists():
            mime = "image/png" if path.suffix.lower() == ".png" else "image/jpeg"
            data = base64.b64encode(path.read_bytes()).decode()
            alt = html.escape(sup.get("photo_alt", "Photo of the creator"))
            photo = f'<img class="avatar" src="data:{mime};base64,{data}" alt="{alt}" width="96" height="96">'
        else:
            print(f"  ! support photo {path} not found; leaving it out", file=sys.stderr)
    title = html.escape(sup.get("title") or "Support this tool")
    if sup.get("linkedin"):
        title += (f' <a class="li" href="{html.escape(sup["linkedin"])}" target="_blank" rel="noopener noreferrer" '
                  f'aria-label="LinkedIn profile" title="LinkedIn">{icon("linkedin")}</a>')
    paras = sup.get("message") or [cfg.get("donate_message", "If this tool helped you, you can support it.")]
    if isinstance(paras, str):
        paras = [paras]
    body = "".join(f"<p>{html.escape(p)}</p>" for p in paras)
    return (f'<section class="donate no-print card">{photo}<div class="about"><h2>{title}</h2>{body}'
            f'<div class="links">{"".join(links)}</div></div></section>')


def features_html(data):
    """The "What it checks" cards, built here rather than in the page's JavaScript so
    search engines and link previews that don't run scripts can read them."""
    cat = {c["key"]: c for c in data["categories"]}
    per_cat, common = {}, 0
    for s in data["snps"]:
        per_cat[s["cat"]] = per_cat.get(s["cat"], 0) + 1
        common += cat[s["cat"]]["group"] == "Common"
    cards = [(f"{per_cat.get(k, 0)} SNPs", title, text) for k, title, text in [
        ("thy_conv", "T4 → T3 conversion", "DIO1 and DIO2, including Thr92Ala."),
        ("thy_tsh", "TSH set-point", "PDE8B, FOXE1 and VAV3."),
        ("thy_auto", "Hashimoto's & Graves'", "PTPN22, CTLA4, TSHR and celiac-linked HLA-DQ."),
        ("thy_cancer", "Nodules & cancer", "FOXE1, PTCSC3 and DIRC3."),
        ("thy_cof", "Thyroid nutrients", "Selenium, iron, vitamin A/D and B12 genes."),
    ]]
    cards.append((f"{common} SNPs", "Common variants", "MTHFR, COMT, clotting, drug response and more."))
    cards.append(("CSV export", "Every thyroid-gene SNP", f"All SNPs your chip has in {len(data['genes'])} thyroid genes."))
    return "".join(f'<div class="card feat"><div class="n">{html.escape(n)}</div><h3>{html.escape(h)}</h3>'
                   f'<p>{html.escape(t)}</p></div>' for n, h, t in cards)


def seo_html(cfg, data, page):
    """Canonical link, social-preview (Open Graph / Twitter) tags, structured data and
    search-console verification tags, all derived from the page's own title/description."""
    url = data["site"]["url"]
    title = html.unescape(re.search(r"<title>(.*?)</title>", page).group(1))
    desc = html.unescape(re.search(r'<meta name="description" content="(.*?)">', page).group(1))
    a = lambda v: html.escape(v, quote=True)
    tags = [f'<link rel="canonical" href="{a(url)}">',
            '<meta name="robots" content="index,follow">',
            '<meta property="og:type" content="website">',
            '<meta property="og:site_name" content="Thyroid SNP Check">',
            f'<meta property="og:title" content="{a(title)}">',
            f'<meta property="og:description" content="{a(desc)}">',
            f'<meta property="og:url" content="{a(url)}">']
    card = "summary"
    if (DOCS / OG_IMAGE).exists():
        img = url.rstrip("/") + "/" + OG_IMAGE
        card = "summary_large_image"
        tags += [f'<meta property="og:image" content="{a(img)}">',
                 '<meta property="og:image:width" content="1200">',
                 '<meta property="og:image:height" content="630">',
                 '<meta property="og:image:alt" content="Thyroid SNP Check: a free thyroid gene report from your DNA raw data">',
                 f'<meta name="twitter:image" content="{a(img)}">']
    tags.insert(7, f'<meta name="twitter:card" content="{card}">')
    for engine, name in (("google", "google-site-verification"), ("bing", "msvalidate.01")):
        code = ((cfg.get("verification") or {}).get(engine) or "").strip()
        if code:
            tags.append(f'<meta name="{name}" content="{a(code)}">')
    ld = {
        "@context": "https://schema.org", "@type": "WebApplication",
        "name": "Thyroid SNP Check", "url": url, "description": desc,
        "applicationCategory": "HealthApplication", "operatingSystem": "Any (runs in a web browser)",
        "browserRequirements": "Requires JavaScript",
        "isAccessibleForFree": True, "offers": {"@type": "Offer", "price": "0", "priceCurrency": "USD"},
        "license": "https://polyformproject.org/licenses/noncommercial/1.0.0/",
        "codeRepository": data["site"]["repo"], "inLanguage": "en",
    }
    tags.append('<script type="application/ld+json">'
                + json.dumps(ld, ensure_ascii=False).replace("</", "<\\/") + "</script>")
    return "\n".join(tags)


def apply_analytics(page, cfg):
    a = cfg.get("analytics") or {}
    provider, ident = (a.get("provider") or "").strip().lower(), (a.get("id") or "").strip()
    if not provider:
        return page.replace("<!--__ANALYTICS__-->", ""), None
    if provider not in ANALYTICS:
        sys.exit(f"Unknown analytics provider {provider!r}; use one of {', '.join(ANALYTICS)}")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", ident):
        sys.exit("analytics.id must be the bare token / site code (letters, digits, - or _)")
    spec = {k: v.format(id=ident) if isinstance(v, str) else v for k, v in ANALYTICS[provider].items()}
    csp_old = "script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data: blob:; connect-src 'none';"
    assert csp_old in page, "CSP in template changed; update build_site.py"
    img = f" {spec['img']}" if spec["img"] else ""
    csp_new = (f"script-src 'unsafe-inline' {spec['script']}; style-src 'unsafe-inline'; "
               f"img-src data: blob:{img}; connect-src {spec['connect']};")
    page = page.replace(csp_old, csp_new)
    page = page.replace("<!--__ANALYTICS__-->", spec["tag"])
    for old, new in PRIVACY_WITH_ANALYTICS:
        assert old in page, f"privacy text changed in template: {old[:50]}"
        page = page.replace(old, new.format(name=spec["name"]))
    return page, spec["name"]


def build(cfg, analytics=True, seo=True):
    coords = json.loads((ROOT / "gene_coords_grch37.json").read_text(encoding="utf-8"))
    data = {
        "categories": [{"key": k, "title": t, "group": g} for k, t, g in CATEGORIES],
        "snps": SNPS,
        "apoe": {"|".join(k): v for k, v in APOE_TABLE.items()},
        "apoeNotes": APOE_NOTES,
        "genes": {g: coords[g] for g in THYROID_GENES if g in coords},
        "positions": load_positions(),
        "site": {"repo": cfg.get("repo_url") or DEFAULT_REPO, "url": cfg.get("site_url") or DEFAULT_SITE},
    }
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    page = (SITE / "template.html").read_text(encoding="utf-8")
    page = page.replace("/*__SNP_DATA__*/", payload).replace("<!--__DONATE__-->", support_html(cfg))
    page = page.replace("<!--__FEATURES__-->", features_html(data))
    page = page.replace("<!--__SEO__-->", seo_html(cfg, data, page) if seo else "")
    used = None
    if analytics:
        page, used = apply_analytics(page, cfg)
    else:
        page = page.replace("<!--__ANALYTICS__-->", "")
    return page, data, used


def to_preview(page):
    """claude.ai Artifact variant: the viewer supplies the document skeleton and its own
    sandbox (which blocks downloads and printing), so drop our skeleton/CSP and the
    download buttons, and say where downloads are available."""
    page = re.sub(r"<!doctype html>\s*<html[^>]*>\s*<head>", "", page, flags=re.I)
    page = re.sub(r'<meta charset="utf-8">\s*<meta name="viewport"[^>]*>', "", page)
    page = re.sub(r"<!-- The page cannot make.*?-->\s*", "", page, flags=re.S)
    page = re.sub(r'<meta http-equiv="Content-Security-Policy".*?>', "", page, flags=re.S)
    page = re.sub(r"</head>\s*<body>", "", page)
    page = re.sub(r"</body>\s*</html>\s*$", "", page)
    page = page.replace(".preview-note{display:none}",
                        ".preview-note{display:block}.dl-only{display:none!important}")
    page = page.replace('<div id="report"></div>',
                        '<p class="preview-note small">Preview build: downloads and PDF export are '
                        'turned off here and work on the hosted site.</p><div id="report"></div>')
    return page


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preview", help="also write an Artifact preview to this path")
    args = ap.parse_args()
    cfg, cfg_name = load_config()
    page, data, used = build(cfg)
    DOCS.mkdir(exist_ok=True)
    (DOCS / ".nojekyll").touch()
    (DOCS / "index.html").write_text(page, encoding="utf-8")
    (DOCS / "sitemap.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
        f'  <url><loc>{html.escape(data["site"]["url"])}</loc></url>\n'
        '</urlset>\n', encoding="utf-8")
    print(f"Wrote {DOCS / 'index.html'} ({len(page) / 1024:.0f} KB, {len(SNPS)} SNPs, "
          f"{len(data['genes'])} genes, {len(data['positions'])} SNP positions; "
          f"config: {cfg_name or 'none'}; analytics: {used or 'off'})")
    if args.preview:
        preview, _, _ = build(cfg, analytics=False, seo=False)
        Path(args.preview).write_text(to_preview(preview), encoding="utf-8")
        print(f"Wrote preview {args.preview}")


if __name__ == "__main__":
    main()
