#!/usr/bin/env python3
"""
Draw docs/og-image.png, the 1200x630 card that social sites and chat apps show when
someone shares the link. Needs Pillow; only rerun it when the wording should change.

    python site/make_og_image.py
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parent.parent / "docs" / "og-image.png"
W, H, S = 1200, 630, 2  # drawn at 2x, then scaled down for smooth edges
BG, INK, MUTED, ACCENT, PILL = "#f7f7f5", "#1d1d1f", "#5b5b60", "#2f5d8a", "#e4edf6"
FONTS = Path("C:/Windows/Fonts")


def font(bold, size):
    for name in (["segoeuib.ttf", "arialbd.ttf"] if bold else ["segoeui.ttf", "arial.ttf"]):
        if (FONTS / name).exists():
            return ImageFont.truetype(str(FONTS / name), size * S)
    return ImageFont.load_default(size * S)


def logo(x, y, scale):
    """The site's thyroid mark (the favicon's SVG path), as a polygon."""
    # M16 6 c-3 0-4 3-6 3 S4 8 4 14 s4 10 8 10 c2 0 3-2 4-4 1 2 2 4 4 4 4 0 8-4 8-10 s-4-5-6-5-3-3-6-3z
    segs = [((16, 6), (13, 6), (12, 9), (10, 9)), ((10, 9), (8, 9), (4, 8), (4, 14)),
            ((4, 14), (4, 20), (8, 24), (12, 24)), ((12, 24), (14, 24), (15, 22), (16, 20)),
            ((16, 20), (17, 22), (18, 24), (20, 24)), ((20, 24), (24, 24), (28, 20), (28, 14)),
            ((28, 14), (28, 8), (24, 9), (22, 9)), ((22, 9), (19, 9), (19, 6), (16, 6))]
    pts = []
    for p0, p1, p2, p3 in segs:
        for i in range(24):
            t = i / 24
            pts.append(tuple((1 - t) ** 3 * a + 3 * (1 - t) ** 2 * t * b + 3 * (1 - t) * t * t * c + t ** 3 * d
                             for a, b, c, d in zip(p0, p1, p2, p3)))
    return [((x + px * scale) * S, (y + py * scale) * S) for px, py in pts]


def main():
    img = Image.new("RGB", (W * S, H * S), BG)
    d = ImageDraw.Draw(img)
    d.rectangle([0, 0, W * S, 14 * S], fill=ACCENT)
    d.polygon(logo(72, 62, 2.4), fill=ACCENT)
    d.text((158 * S, 70 * S), "Thyroid SNP Check", font=font(True, 38), fill=ACCENT)
    y = 178
    for line in ("What does your DNA say", "about your thyroid?"):
        d.text((72 * S, y * S), line, font=font(True, 72), fill=INK)
        y += 88
    d.text((72 * S, (y + 22) * S), "Free thyroid gene report from 23andMe, AncestryDNA,",
           font=font(False, 34), fill=MUTED)
    d.text((72 * S, (y + 68) * S), "MyHeritage and FamilyTreeDNA raw data.", font=font(False, 34), fill=MUTED)
    pill = "Runs on your device  ·  nothing is uploaded"
    f = font(True, 28)
    w = d.textlength(pill, font=f) / S
    d.rounded_rectangle([72 * S, 520 * S, (72 + w + 48) * S, 572 * S], radius=26 * S, fill=PILL)
    d.text(((72 + 24) * S, 546 * S), pill, font=f, fill=ACCENT, anchor="lm")
    OUT.parent.mkdir(exist_ok=True)
    img.resize((W, H), Image.LANCZOS).save(OUT, optimize=True)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
