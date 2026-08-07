#!/usr/bin/env python3
"""
Assemble main Figure 2 (2x2 panel layout) for publication.

Panels
------
a  Overview / summary figure
b  ABO hotspot
c  Chromosome / novel QTL overview
d  Novel pQTL MAF scatter

Dependencies
------------
    pip install matplotlib numpy pillow pymupdf

Usage
-----
    # Default: look for panel images under ./input relative to this script
    python plot_Fig2_main.py

    # Or specify panel files and output directory
    python plot_Fig2_main.py \\
        --panel-a input/Fig2A.jpg \\
        --panel-b input/hotspot_ABO_2.pdf \\
        --panel-c input/Fig2C_update.png \\
        --panel-d input/novel_pQTL_MAF_scatter_final.png \\
        --outdir output

Outputs
-------
    Fig2_main_2x2.png
    Fig2_main_2x2.pdf   (scaled to fit one A4 page)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def load_image(path: Path) -> np.ndarray:
    """Load jpg/png directly; render the first page for PDF."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Panel image not found: {path}")

    suffix = path.suffix.lower()
    if suffix == ".pdf":
        import fitz  # PyMuPDF

        doc = fitz.open(path)
        page = doc[0]
        pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
        img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        doc.close()
    else:
        img = Image.open(path).convert("RGB")

    return np.asarray(img)


def assemble_figure2(
    panel_paths: list[Path],
    output_png: Path,
    output_pdf: Path,
    scales: list[float] | None = None,
    figsize: tuple[float, float] = (16, 12),
    dpi: int = 300,
) -> None:
    """
    Compose four panels into a 2x2 figure.

    Layout
    ------
    Top row width ratio a:b = 6:4
    Bottom row width ratio c:d = 5:5
    """
    if len(panel_paths) != 4:
        raise ValueError("Exactly four panel paths (a–d) are required.")

    if scales is None:
        # Panel b slightly enlarged for visual balance
        scales = [1.0, 1.18, 1.0, 1.0]
    if len(scales) != 4:
        raise ValueError("scales must contain four values.")

    labels = ["a", "b", "c", "d"]
    images = [load_image(p) for p in panel_paths]

    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = fig.add_gridspec(2, 10)
    axes = [
        fig.add_subplot(gs[0, :6]),
        fig.add_subplot(gs[0, 6:]),
        fig.add_subplot(gs[1, :5]),
        fig.add_subplot(gs[1, 5:]),
    ]

    label_x = -0.08
    label_y = 1.08

    for i, (ax, img, label, scale) in enumerate(zip(axes, images, labels, scales)):
        h, w = img.shape[:2]
        ax.imshow(img)

        # Keep aspect ratio; allow per-panel zoom
        cx, cy = w / 2, h / 2
        half_w = (w / 2) / scale
        half_h = (h / 2) / scale
        ax.set_xlim(cx - half_w, cx + half_w)
        ax.set_ylim(cy + half_h, cy - half_h)
        ax.axis("off")

        # Panel a label is placed in figure coordinates after layout
        if i == 0:
            continue

        ax.text(
            label_x,
            label_y,
            label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=22,
            fontweight="bold",
            color="black",
            clip_on=False,
        )

    plt.subplots_adjust(
        left=0.04, right=0.98, top=0.96, bottom=0.03, wspace=0.05, hspace=0.10
    )

    # Align panel-a label with c (x) and b (y)
    pos_b = axes[1].get_position()
    pos_c = axes[2].get_position()
    fig_x_a = pos_c.x0 + label_x * pos_c.width
    fig_y_a = pos_b.y0 + label_y * pos_b.height
    fig.text(
        fig_x_a,
        fig_y_a,
        "a",
        transform=fig.transFigure,
        ha="left",
        va="top",
        fontsize=22,
        fontweight="bold",
        color="black",
        clip_on=False,
    )

    output_png = Path(output_png)
    output_pdf = Path(output_pdf)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_png, dpi=dpi)

    # Scale PDF to fit within one A4 page without changing layout proportions
    a4_w_in, a4_h_in = 210 / 25.4, 297 / 25.4
    fig_w, fig_h = fig.get_size_inches()
    fit_scale = min(a4_w_in / fig_w, a4_h_in / fig_h)
    fig.set_size_inches(fig_w * fit_scale, fig_h * fit_scale)
    fig.savefig(output_pdf, dpi=dpi, format="pdf", pad_inches=0)
    fig.set_size_inches(fig_w, fig_h)
    plt.close(fig)

    print(f"Saved: {output_png}")
    print(f"Saved: {output_pdf}")
    print(
        f"PDF page: {fig_w * fit_scale * 25.4:.1f} x {fig_h * fit_scale * 25.4:.1f} mm "
        f"(fits in A4 210 x 297 mm)"
    )


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    default_input = here / "input"

    parser = argparse.ArgumentParser(
        description="Assemble publication main Figure 2 (2x2 panels)."
    )
    parser.add_argument(
        "--panel-a",
        type=Path,
        default=default_input / "Fig2A.jpg",
        help="Path to panel a image (jpg/png/pdf).",
    )
    parser.add_argument(
        "--panel-b",
        type=Path,
        default=default_input / "hotspot_ABO_2.pdf",
        help="Path to panel b image (jpg/png/pdf).",
    )
    parser.add_argument(
        "--panel-c",
        type=Path,
        default=default_input / "Fig2C_update.png",
        help="Path to panel c image (jpg/png/pdf).",
    )
    parser.add_argument(
        "--panel-d",
        type=Path,
        default=default_input / "novel_pQTL_MAF_scatter_final.png",
        help="Path to panel d image (jpg/png/pdf).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=here / "output",
        help="Directory for Fig2_main_2x2.png/pdf.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    assemble_figure2(
        panel_paths=[args.panel_a, args.panel_b, args.panel_c, args.panel_d],
        output_png=outdir / "Fig2_main_2x2.png",
        output_pdf=outdir / "Fig2_main_2x2.pdf",
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
