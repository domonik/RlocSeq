#!/usr/bin/env python
"""Plot the WT and delta-rbp3 psaA FISH localisation result.

The supplied ``distributions.tsv`` already contains the radial intensity
measurements.  As in the repository's FISH plotting scripts, each spot profile
is min-max normalized independently and the normalized profiles are averaged
at every radial distance.  This figure intentionally uses all 30 supplied
rings (0--754 nm) to reproduce the attached reference; the current general
pipeline figures truncate to rings 0--20 before normalization.

The final SVG embeds the prepared representative TIFF panels and draws exact,
condition-specific vector intensity bars.  Inkscape is used for the
publication-resolution PNG export, matching the existing Snakemake workflow.
"""

from __future__ import annotations

import argparse
import base64
import io
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from PIL import Image
import tifffile


SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)

RING_COUNT = 30
NM_PER_RING = 26
GROUPS = ("WT", "rbp3_C2")
DISPLAY_LABELS = {"WT": "WT", "rbp3_C2": "Δrbp3"}
COLORS = {"WT": "#689BF7", "rbp3_C2": "#4671F2"}
LINE_DASH = {"WT": None, "rbp3_C2": "10 6"}
COLORBAR_MAX = {"WT": 1800, "rbp3_C2": 1500}
REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ProfileSummary:
    treatment: str
    display_label: str
    distances_nm: np.ndarray
    mean_profile: np.ndarray
    n_spots: int
    n_images: int
    n_cells: int
    peak_distance_nm: int


def _svg_tag(name: str) -> str:
    return f"{{{SVG_NS}}}{name}"


def _parse_args() -> argparse.Namespace:
    root = REPO_ROOT / "initialData/Rbp3_psaA_Localisation"
    output_dir = REPO_ROOT / "Pipeline/FISH/Rbp3PsaALocalisation"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", type=Path, default=root / "distributions.tsv")
    parser.add_argument(
        "--wt-image",
        type=Path,
        default=root / "Composite_Image 9_Out-Channel Alignment-04_4.916-1_WT.tif",
    )
    parser.add_argument(
        "--delta-image",
        type=Path,
        default=root / "Composite_Image 18_Out-Channel Alignment-03_5.323-1_drbp3.tif",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=output_dir / "Rbp3PsaALocalisation.svg",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=output_dir / "Rbp3PsaALocalisation.png",
    )
    parser.add_argument(
        "--output-data",
        type=Path,
        default=output_dir / "PlotData.tsv",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG export resolution (default: 300 dpi)",
    )
    return parser.parse_args()


def _require_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input file(s): " + ", ".join(missing))


def _summarize_profiles(distribution_file: Path) -> dict[str, ProfileSummary]:
    df = pd.read_csv(distribution_file, sep="\t")
    ring_columns = [f"ring_{idx}_mean" for idx in range(RING_COUNT)]
    required = {"Sonde", "Treatment", "Image", "Cell", *ring_columns}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Distribution table lacks required column(s): {missing}")

    df = df.loc[df["Sonde"].eq("psaA") & df["Treatment"].isin(GROUPS)].copy()
    summaries: dict[str, ProfileSummary] = {}
    distances = np.arange(RING_COUNT, dtype=int) * NM_PER_RING

    for treatment in GROUPS:
        group = df.loc[df["Treatment"].eq(treatment)]
        if group.empty:
            raise ValueError(f"No psaA profiles found for Treatment={treatment!r}")

        values = group.loc[:, ring_columns].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            bad_count = int((~np.isfinite(values)).sum())
            raise ValueError(
                f"Treatment={treatment!r} contains {bad_count} non-finite ring values"
            )

        row_min = values.min(axis=1, keepdims=True)
        row_span = values.max(axis=1, keepdims=True) - row_min
        zero_span = np.flatnonzero(row_span[:, 0] <= 0)
        if zero_span.size:
            raise ValueError(
                f"Treatment={treatment!r} contains {zero_span.size} constant profiles"
            )

        normalized = (values - row_min) / row_span
        mean_profile = normalized.mean(axis=0)
        peak_index = int(np.argmax(mean_profile))
        summaries[treatment] = ProfileSummary(
            treatment=treatment,
            display_label=DISPLAY_LABELS[treatment],
            distances_nm=distances.copy(),
            mean_profile=mean_profile,
            n_spots=len(group),
            n_images=int(group["Image"].nunique()),
            n_cells=int(group["Cell"].nunique()),
            peak_distance_nm=int(distances[peak_index]),
        )

    return summaries


def _write_plot_data(
    summaries: dict[str, ProfileSummary], output_file: Path
) -> None:
    frames = []
    normalization = (
        "per-spot min-max across all 30 supplied rings "
        "(ring_0_mean..ring_29_mean)"
    )
    for treatment in GROUPS:
        summary = summaries[treatment]
        frames.append(
            pd.DataFrame(
                {
                    "Probe": "psaA",
                    "Treatment": summary.treatment,
                    "DisplayLabel": summary.display_label,
                    "Ring": np.arange(RING_COUNT, dtype=int),
                    "Distance_nm": summary.distances_nm,
                    "MeanNormalizedIntensity": summary.mean_profile,
                    "n_spots": summary.n_spots,
                    "n_images": summary.n_images,
                    "n_cells": summary.n_cells,
                    "PeakDistance_nm": summary.peak_distance_nm,
                    "Normalization": normalization,
                }
            )
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(
        output_file, sep="\t", index=False, float_format="%.15g"
    )


def _png_data_uri(tiff_file: Path) -> tuple[str, tuple[int, int]]:
    array = tifffile.imread(tiff_file)
    if array.dtype != np.uint8:
        raise ValueError(f"Expected uint8 image data in {tiff_file}, got {array.dtype}")
    if array.ndim != 3 or array.shape[-1] not in (3, 4):
        raise ValueError(
            f"Expected an RGB/RGBA image in {tiff_file}, got shape {array.shape}"
        )

    image = Image.fromarray(array)
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}", (image.width, image.height)


def _add_text(
    parent: ET.Element,
    text: str,
    x: float,
    y: float,
    **attributes: str,
) -> ET.Element:
    attrs = {"x": f"{x:g}", "y": f"{y:g}", **attributes}
    node = ET.SubElement(parent, _svg_tag("text"), attrs)
    node.text = text
    return node


def _write_svg(
    summaries: dict[str, ProfileSummary],
    wt_image: Path,
    delta_image: Path,
    output_file: Path,
) -> None:
    width, height = 1200, 400
    image_y, image_size = 45, 286
    wt_image_x, delta_image_x = 20, 408
    wt_bar_x, delta_bar_x = 314, 702
    bar_y, bar_height = 114, 148
    plot_left, plot_right = 830, 1170
    plot_top, plot_bottom = 45, 331
    plot_width = plot_right - plot_left
    plot_height = plot_bottom - plot_top
    x_max = 780.0
    y_max = 1.0

    root = ET.Element(
        _svg_tag("svg"),
        {
            "version": "1.1",
            "width": str(width),
            "height": str(height),
            "viewBox": f"0 0 {width} {height}",
            "role": "img",
            "aria-labelledby": "figure-title figure-description",
        },
    )
    title = ET.SubElement(root, _svg_tag("title"), {"id": "figure-title"})
    title.text = "psaA FISH localisation in WT and delta-rbp3"
    description = ET.SubElement(
        root, _svg_tag("desc"), {"id": "figure-description"}
    )
    description.text = (
        "Representative WT and delta-rbp3 FISH images with separate intensity "
        "color bars and a combined mean normalized radial intensity plot. "
        "Both white image scale bars represent 2 micrometres."
    )

    defs = ET.SubElement(root, _svg_tag("defs"))
    style = ET.SubElement(defs, _svg_tag("style"))
    style.text = (
        "text { font-family: Arial, Helvetica, sans-serif; fill: #000; } "
        ".axis { stroke: #000; stroke-width: 1.5; shape-rendering: geometricPrecision; } "
        ".tick { stroke: #000; stroke-width: 1.25; }"
    )
    clip = ET.SubElement(defs, _svg_tag("clipPath"), {"id": "plot-clip"})
    ET.SubElement(
        clip,
        _svg_tag("rect"),
        {
            "x": str(plot_left),
            "y": str(plot_top),
            "width": str(plot_width),
            "height": str(plot_height),
        },
    )
    gradient = ET.SubElement(
        defs,
        _svg_tag("linearGradient"),
        {"id": "black-green", "x1": "0", "y1": "1", "x2": "0", "y2": "0"},
    )
    ET.SubElement(
        gradient,
        _svg_tag("stop"),
        {"offset": "0", "stop-color": "#000000"},
    )
    ET.SubElement(
        gradient,
        _svg_tag("stop"),
        {"offset": "1", "stop-color": "#00ff00"},
    )

    ET.SubElement(
        root,
        _svg_tag("rect"),
        {"x": "0", "y": "0", "width": str(width), "height": str(height), "fill": "#fff"},
    )

    panel_specs = (
        ("WT", wt_image, wt_image_x, wt_bar_x),
        ("rbp3_C2", delta_image, delta_image_x, delta_bar_x),
    )
    for treatment, image_file, image_x, colorbar_x in panel_specs:
        image_uri, image_dimensions = _png_data_uri(image_file)
        if image_dimensions != (286, 286):
            raise ValueError(
                f"Expected a 286 x 286 representative image, got {image_dimensions}: "
                f"{image_file}"
            )
        ET.SubElement(
            root,
            _svg_tag("image"),
            {
                "x": str(image_x),
                "y": str(image_y),
                "width": str(image_size),
                "height": str(image_size),
                "href": image_uri,
                "preserveAspectRatio": "xMidYMid meet",
            },
        )
        bar_rect_width = 18
        ET.SubElement(
            root,
            _svg_tag("rect"),
            {
                "x": str(colorbar_x),
                "y": str(bar_y),
                "width": str(bar_rect_width),
                "height": str(bar_height),
                "fill": "url(#black-green)",
                "stroke": "#000",
                "stroke-width": "1.25",
            },
        )
        maximum = COLORBAR_MAX[treatment]
        for tick_value in np.linspace(0, maximum, 5, dtype=int):
            tick_y = bar_y + bar_height * (1 - tick_value / maximum)
            ET.SubElement(
                root,
                _svg_tag("line"),
                {
                    "x1": str(colorbar_x + bar_rect_width),
                    "x2": str(colorbar_x + bar_rect_width + 5),
                    "y1": f"{tick_y:.3f}",
                    "y2": f"{tick_y:.3f}",
                    "class": "tick",
                },
            )
            _add_text(
                root,
                str(tick_value),
                colorbar_x + bar_rect_width + 10,
                tick_y + 5,
                **{"font-size": "14"},
            )
        _add_text(
            root,
            DISPLAY_LABELS[treatment],
            image_x + image_size / 2,
            27,
            **{
                "text-anchor": "middle",
                "font-size": "18",
                "font-style": "italic" if treatment == "rbp3_C2" else "normal",
            },
        )

    _add_text(
        root,
        "psaA",
        (plot_left + plot_right) / 2,
        27,
        **{"text-anchor": "middle", "font-size": "18", "font-style": "italic"},
    )

    ET.SubElement(
        root,
        _svg_tag("line"),
        {
            "x1": str(plot_left),
            "x2": str(plot_left),
            "y1": str(plot_top),
            "y2": str(plot_bottom),
            "class": "axis",
        },
    )
    ET.SubElement(
        root,
        _svg_tag("line"),
        {
            "x1": str(plot_left),
            "x2": str(plot_right),
            "y1": str(plot_bottom),
            "y2": str(plot_bottom),
            "class": "axis",
        },
    )

    def x_to_svg(value: float) -> float:
        return plot_left + (value / x_max) * plot_width

    def y_to_svg(value: float) -> float:
        return plot_bottom - (value / y_max) * plot_height

    for value in (0, 250, 500, 750):
        x = x_to_svg(value)
        ET.SubElement(
            root,
            _svg_tag("line"),
            {
                "x1": f"{x:.3f}",
                "x2": f"{x:.3f}",
                "y1": str(plot_bottom),
                "y2": str(plot_bottom + 6),
                "class": "tick",
            },
        )
        _add_text(
            root,
            str(value),
            x,
            plot_bottom + 23,
            **{"text-anchor": "middle", "font-size": "14"},
        )

    for value in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = y_to_svg(value)
        ET.SubElement(
            root,
            _svg_tag("line"),
            {
                "x1": str(plot_left - 6),
                "x2": str(plot_left),
                "y1": f"{y:.3f}",
                "y2": f"{y:.3f}",
                "class": "tick",
            },
        )
        _add_text(
            root,
            f"{value:.2f}",
            plot_left - 11,
            y + 5,
            **{"text-anchor": "end", "font-size": "14"},
        )

    _add_text(
        root,
        "Distance to probe [nm]",
        (plot_left + plot_right) / 2,
        384,
        **{"text-anchor": "middle", "font-size": "16", "font-weight": "bold"},
    )
    _add_text(
        root,
        "Norm. intensity",
        772,
        (plot_top + plot_bottom) / 2,
        **{
            "text-anchor": "middle",
            "font-size": "16",
            "font-weight": "bold",
            "transform": f"rotate(-90 772 {(plot_top + plot_bottom) / 2:g})",
        },
    )

    plot_group = ET.SubElement(
        root, _svg_tag("g"), {"clip-path": "url(#plot-clip)"}
    )
    for treatment in GROUPS:
        summary = summaries[treatment]
        peak_x = x_to_svg(summary.peak_distance_nm)
        ET.SubElement(
            plot_group,
            _svg_tag("line"),
            {
                "x1": f"{peak_x:.3f}",
                "x2": f"{peak_x:.3f}",
                "y1": str(plot_top),
                "y2": str(plot_bottom),
                "stroke": COLORS[treatment],
                "stroke-width": "1.5",
                "stroke-dasharray": "5 5",
                "opacity": "0.9",
            },
        )

    for treatment in GROUPS:
        summary = summaries[treatment]
        points = " ".join(
            f"{x_to_svg(distance):.3f},{y_to_svg(value):.3f}"
            for distance, value in zip(
                summary.distances_nm, summary.mean_profile, strict=True
            )
        )
        ET.SubElement(
            plot_group,
            _svg_tag("polyline"),
            {
                "points": points,
                "fill": "none",
                "stroke": COLORS[treatment],
                "stroke-width": "3",
                "stroke-linecap": "round",
                "stroke-linejoin": "round",
                **(
                    {"stroke-dasharray": LINE_DASH[treatment]}
                    if LINE_DASH[treatment]
                    else {}
                ),
            },
        )
        peak_index = int(np.argmax(summary.mean_profile))
        ET.SubElement(
            plot_group,
            _svg_tag("circle"),
            {
                "cx": f"{x_to_svg(summary.peak_distance_nm):.3f}",
                "cy": f"{y_to_svg(summary.mean_profile[peak_index]):.3f}",
                "r": "3.5",
                "fill": COLORS[treatment],
            },
        )

    peak_label_settings = {
        "WT": {"dx": -4, "y": plot_bottom - 7, "anchor": "end"},
        "rbp3_C2": {"dx": 4, "y": plot_bottom - 25, "anchor": "start"},
    }
    for treatment in GROUPS:
        summary = summaries[treatment]
        settings = peak_label_settings[treatment]
        _add_text(
            root,
            str(summary.peak_distance_nm),
            x_to_svg(summary.peak_distance_nm) + settings["dx"],
            settings["y"],
            **{
                "text-anchor": settings["anchor"],
                "font-size": "12",
                "font-weight": "bold",
                "fill": COLORS[treatment],
                "style": (
                    f"fill:{COLORS[treatment]};paint-order:stroke;stroke:#fff;"
                    "stroke-width:3px;stroke-linejoin:round"
                ),
            },
        )

    legend_x, legend_y = plot_right - 158, plot_top + 10
    ET.SubElement(
        root,
        _svg_tag("rect"),
        {
            "x": str(legend_x - 8),
            "y": str(legend_y - 7),
            "width": "166",
            "height": "55",
            "rx": "3",
            "fill": "#fff",
            "opacity": "0.9",
        },
    )
    for row, treatment in enumerate(GROUPS):
        summary = summaries[treatment]
        y = legend_y + row * 25
        ET.SubElement(
            root,
            _svg_tag("line"),
            {
                "x1": str(legend_x),
                "x2": str(legend_x + 30),
                "y1": str(y),
                "y2": str(y),
                "stroke": COLORS[treatment],
                "stroke-width": "3",
                "stroke-linecap": "round",
                **(
                    {"stroke-dasharray": LINE_DASH[treatment]}
                    if LINE_DASH[treatment]
                    else {}
                ),
            },
        )
        if treatment == "WT":
            _add_text(
                root,
                f"WT, n = {summary.n_spots}",
                legend_x + 39,
                y + 5,
                **{"font-size": "14"},
            )
        else:
            legend_text = ET.SubElement(
                root,
                _svg_tag("text"),
                {
                    "x": str(legend_x + 39),
                    "y": str(y + 5),
                    "font-size": "14",
                },
            )
            gene = ET.SubElement(
                legend_text, _svg_tag("tspan"), {"font-style": "italic"}
            )
            gene.text = "Δrbp3"
            suffix = ET.SubElement(legend_text, _svg_tag("tspan"))
            suffix.text = f"\u00a0#2, n = {summary.n_spots}"

    output_file.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(output_file, encoding="utf-8", xml_declaration=True)


def _export_png(svg_file: Path, png_file: Path, dpi: int) -> None:
    if dpi <= 0:
        raise ValueError(f"DPI must be positive, got {dpi}")
    inkscape = shutil.which("inkscape")
    if inkscape is None:
        raise RuntimeError("Inkscape is required to export the PNG")
    png_file.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            inkscape,
            str(svg_file),
            "--export-type=png",
            f"--export-filename={png_file}",
            f"--export-dpi={dpi}",
        ],
        check=True,
    )


def main() -> None:
    args = _parse_args()
    _require_files(
        (
            args.input_tsv,
            args.wt_image,
            args.delta_image,
        )
    )
    summaries = _summarize_profiles(args.input_tsv)
    _write_plot_data(summaries, args.output_data)
    _write_svg(
        summaries,
        args.wt_image,
        args.delta_image,
        args.output_svg,
    )
    _export_png(args.output_svg, args.output_png, args.dpi)

    for treatment in GROUPS:
        summary = summaries[treatment]
        print(
            f"{summary.display_label}: n={summary.n_spots}, "
            f"images={summary.n_images}, cells={summary.n_cells}, "
            f"peak={summary.peak_distance_nm} nm"
        )
    print(f"Wrote {args.output_svg}")
    print(f"Wrote {args.output_png}")
    print(f"Wrote {args.output_data}")


if __name__ == "__main__":
    main()
