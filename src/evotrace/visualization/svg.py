"""Dependency-free, scalable conservation profile."""

from pathlib import Path


def profile_svg(columns, out: Path, key="conservation_score", title="Conservation profile"):
    w, h, pad = 1000, 340, 55
    vals = [float(x[key]) for x in columns]
    coords = []
    for i, v in enumerate(vals):
        x = pad + (w - 2 * pad) * (i / max(1, len(vals) - 1))
        y = h - pad - v * (h - 2 * pad)
        coords.append("{:.1f},{:.1f}".format(x, y))
    poly = " ".join(coords)
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100%" viewBox="0 0 {} {}"><rect width="100%" height="100%" fill="white"/><text x="{}" y="28" font-family="sans-serif" font-size="18">{}</text><line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#777"/><polyline fill="none" stroke="#147d8a" stroke-width="2" points="{}"/><text x="{}" y="{}" font-family="sans-serif" font-size="12">Position</text></svg>'.format(
        w, h, w / 2, title, pad, h - pad, w - pad, h - pad, poly, w / 2 - 20, h - 8
    )
    out.write_text(svg, encoding="utf-8")
