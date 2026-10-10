# Thyroid SNP Check

Check consumer DNA raw data for thyroid-related and common SNPs, with plain-language
explanations of each result. Free for non-commercial use.

**Live site:** https://aimoony.github.io/thyroid-snp-check/

Supported files: 23andMe (v3–v5), AncestryDNA, MyHeritage, FamilyTreeDNA (Family Finder,
including older and concatenated downloads), Living DNA, Genes for Good, tellmeGen and
Sano, as .txt, .csv, .zip or .gz. Sequencing VCF files (Nebula, Dante) aren't supported yet.

**Privacy:** the website runs entirely in the visitor's browser. The raw-data file is
never uploaded, there is no server, and the page's Content-Security-Policy
(`connect-src 'none'`) blocks it from making any network request.

## What it checks

- **Thyroid:** T4→T3 conversion (DIO1, DIO2), TSH set-point (PDE8B, FOXE1, VAV3),
  autoimmunity (PTPN22, CTLA4, SH2B3, TSHR, HLA-DQ2/DQ8), nodule/cancer markers
  (FOXE1, PTCSC3, DIRC3), and cofactor genes (selenium, iron, vitamin A/D, B12).
- **Common variants:** methylation (MTHFR, COMT…), cardiovascular and clotting,
  metabolic, pharmacogenomics, detox and nutrition/traits.
- **Every SNP** the chip measured in 44 thyroid-pathway genes, as a CSV export.

Each SNP lists its effect allele (on the GRCh37 plus strand), what it's associated
with, and an evidence level. The APOE result is hidden on the website until the
visitor chooses to reveal it.

## Repository layout

| Path | What it is |
| --- | --- |
| `snp_database.py` | The curated SNP list and interpretations, shared by both tools |
| `snp_report.py` | Command-line version: builds an HTML report and CSVs locally |
| `gene_coords_grch37.json` | Thyroid gene coordinates (GRCh37, from Ensembl) |
| `site/template.html` | Website source |
| `site/build_site.py` | Builds `site/index.html` from the template and database |
| `site/config.example.json` | Template for `site/config.json` (donation links, analytics) |
| `snp_positions_grch37.json` | GRCh37 positions of report SNPs, for files that use their own IDs |
| `docs/index.html` | Built, self-contained website, served by GitHub Pages |

## Use it locally

```bash
python snp_report.py path/to/your_raw_data.zip
```

Writes `snp_report.html`, `snp_extracted.csv` and `thyroid_gene_region_snps.csv`
next to the script. Requires Python 3.9+, no extra packages.

## Build and deploy the website

```bash
python site/build_site.py
```

GitHub Pages serves `docs/` from `main`, so committing the rebuilt `docs/index.html`
publishes it. The same file works on any other static host. Rebuild after editing `snp_database.py`, `site/template.html` or
`site/config.json`.

### Configuration

Copy `site/config.example.json` to `site/config.json` (git-ignored) and fill in:

- `donate`: Venmo (`@name`), Cash App (`$tag`), PayPal (`name`) or full `https://`
  links. Empty values are hidden.
- `analytics`: `{"provider": "cloudflare", "id": "<beacon token>"}` or
  `{"provider": "goatcounter", "id": "<site code>"}`. Both are cookie-free. The build
  opens the Content-Security-Policy for that one service only and updates the
  privacy text on the page to say an anonymous page-view count is sent.

Values in the config end up in the public page source, so only put public values there.

## Adding a SNP

Add an entry to `SNPS` in `snp_database.py`. Alleles must be on the GRCh37 plus
strand, which is how 23andMe reports genotypes. Include the source study in the
description, then rebuild the site.

## License

[PolyForm Noncommercial 1.0.0](LICENSE.md). You're free to use, modify and share this
for personal, research, educational and other non-commercial purposes. It may not be
sold or used commercially.

## Disclaimer

For education and research only. This is not medical advice or a diagnostic test.
Consumer genotyping arrays are research-grade; confirm anything clinically important
with a clinical-grade test.
