#!/usr/bin/env python
"""Plot nrtB FISH images and radial intensity profiles.

The figure is arranged by replicate.  Each row contains an NH4 image and
colorbar, a NO3 image and colorbar, and a radial-profile comparison plot.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import tifffile
from scipy import stats

try:
    import pytesseract
except ModuleNotFoundError:  # pragma: no cover - available in the pipeline env
    pytesseract = None

try:
    from scripts.plotlyTrash import replace_square_images
    from scripts.plotly_template import PT7, register_template
except ModuleNotFoundError:  # pragma: no cover - supports direct execution
    from plotlyTrash import replace_square_images
    from plotly_template import PT7, register_template


register_template()
pio.renderers.default = "browser"

PROBE = "nrtB"
CONDITIONS = ("NH4", "NO3")
CONDITION_COLORS = {
    "NH4": "rgb(0, 106, 56)",
    "NO3": "rgb(0, 48, 135)",
}
RING_COLUMNS = [f"ring_{i}_mean" for i in range(21)]
RING_DISTANCE_NM = 26


def _replicate_number(name: str) -> int:
    match = re.search(r"(\d+)$", name)
    if match is None:
        raise ValueError(f"Cannot sort replicate directory: {name}")
    return int(match.group(1))


def _find_representative_image(replicate_dir: Path, condition: str) -> Path:
    matches = sorted(replicate_dir.glob(f"Composite_Image *_{condition}.tif"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one WT {condition} representative image in {replicate_dir}; "
            f"found {len(matches)}"
        )
    return matches[0]


def _read_colorbar_limits(colorbar_path: Path) -> tuple[int, int]:
    image = tifffile.imread(colorbar_path)
    max_value = 1
    if pytesseract is not None:
        text = pytesseract.image_to_string(image)
        values = [int(value) for value in re.findall(r"\b\d+\b", text)]
        if values:
            max_value = max(values)

    nonwhite = np.any(image < 245, axis=-1)
    green_values = image[..., 1][nonwhite]
    max_green = int(green_values.max()) if green_values.size else 255
    return max_value, max_green


def _summarize_profile(group: pd.DataFrame) -> dict[str, object]:
    values = group[RING_COLUMNS].to_numpy(dtype=float)
    row_min = np.nanmin(values, axis=1)
    row_max = np.nanmax(values, axis=1)
    denominator = row_max - row_min
    normalized = np.full_like(values, np.nan, dtype=float)
    np.divide(
        values - row_min[:, None],
        denominator[:, None],
        out=normalized,
        where=denominator[:, None] != 0,
    )

    n = int(len(group))
    mean = np.nanmean(normalized, axis=0)
    std = np.nanstd(normalized, axis=0, ddof=1) if n > 1 else np.zeros(len(RING_COLUMNS))
    sem = std / np.sqrt(n) if n else np.zeros(len(RING_COLUMNS))
    ci95 = sem * stats.t.ppf(0.975, df=n - 1) if n > 1 else np.zeros(len(RING_COLUMNS))
    return {"mean": mean, "ci95": ci95, "n": n}


def _add_colorbar(fig: go.Figure, row: int, col: int, colorbar_path: Path) -> None:
    max_value, max_green = _read_colorbar_limits(colorbar_path)
    y = np.linspace(0, max_value, 256)
    z = np.linspace(0, 1, 256).reshape(-1, 1)
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=[0],
            y=y,
            colorscale=[
                [0, "rgb(0, 0, 0)"],
                [1, f"rgb(0, {max_green}, 0)"],
            ],
            showscale=False,
            hoverinfo="skip",
        ),
        row=row,
        col=col,
    )
    fig.update_xaxes(
        visible=False,
        showgrid=False,
        showline=False,
        zeroline=False,
        row=row,
        col=col,
    )
    fig.update_yaxes(
        range=[0, max_value],
        tickvals=[0, max_value / 2, max_value],
        ticktext=["0", f"{max_value / 2:g}", f"{max_value:g}"],
        showgrid=False,
        showline=False,
        ticks="outside",
        ticklen=2,
        tickfont=dict(size=PT7),
        side="right",
        row=row,
        col=col,
    )


def _axis_domain(fig: go.Figure, axis: str, index: int) -> tuple[float, float]:
    name = axis if index == 1 else f"{axis}{index}"
    return tuple(getattr(fig.layout, name).domain)


def build_figure(
    distribution_file: str | Path,
    image_root: str | Path = f"initialData/AllFISHData/{PROBE}",
) -> tuple[go.Figure, list[Path]]:
    image_root = Path(image_root)
    distribution_file = Path(distribution_file)
    df = pd.read_csv(distribution_file, sep="\t")
    required_columns = {"Mask", "Image", "Probe", "Figure", "Treatment", *RING_COLUMNS}
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {distribution_file}: {sorted(missing)}")

    df = df[
        (df["Probe"] == PROBE)
        & (df["Mask"] == "Measured")
        & df["Treatment"].isin(CONDITIONS)
    ].copy()
    if df.empty:
        raise ValueError(f"No measured {PROBE} data found in {distribution_file}")

    replicates = sorted(df["Figure"].unique(), key=_replicate_number)
    if not replicates:
        raise ValueError("No nrtB replicates found in the distribution table")

    replicate_dirs = {path.name: path for path in image_root.glob("Replicate *") if path.is_dir()}
    missing_replicates = set(replicates).difference(replicate_dirs)
    if missing_replicates:
        raise FileNotFoundError(f"Missing image directories for replicates: {sorted(missing_replicates)}")

    specs = [
        [
            {"type": "image"},
            {"type": "xy"},
            {"type": "image"},
            {"type": "xy"},
            {"type": "xy"},
        ]
        for _ in replicates
    ]
    fig = make_subplots(
        rows=len(replicates),
        cols=5,
        specs=specs,
        row_heights=[1 / len(replicates)] * len(replicates),
        column_widths=[0.31, 0.025, 0.31, 0.025, 0.33],
        horizontal_spacing=0.03,
        vertical_spacing=0.12,
    )
    fig.update_layout(template="simple_white_custom")

    # Keep each colorbar close to its image while leaving room on the
    # colorbar's right for its tick labels and the following panel.
    desired_image_colorbar_gap = 0.004
    extra_colorbar_shift = 0.02
    for row in range(1, len(replicates) + 1):
        for image_col, colorbar_col in ((1, 2), (3, 4)):
            image_axis = (row - 1) * 5 + image_col
            colorbar_axis = (row - 1) * 5 + colorbar_col
            image_domain = _axis_domain(fig, "xaxis", image_axis)
            colorbar_domain = _axis_domain(fig, "xaxis", colorbar_axis)
            current_gap = colorbar_domain[0] - image_domain[1]
            shift = current_gap - desired_image_colorbar_gap + extra_colorbar_shift
            fig.update_xaxes(
                domain=[colorbar_domain[0] - shift, colorbar_domain[1] - shift],
                row=row,
                col=colorbar_col,
            )
        profile_axis = (row - 1) * 5 + 5
        profile_domain = _axis_domain(fig, "xaxis", profile_axis)
        fig.update_xaxes(
            domain=[profile_domain[0] + 0.06, profile_domain[1]],
            row=row,
            col=5,
        )

    images_to_replace: list[Path] = []
    counts: dict[str, list[int]] = {condition: [] for condition in CONDITIONS}

    for row, replicate in enumerate(replicates, start=1):
        replicate_dir = replicate_dirs[replicate]
        for condition, image_col, colorbar_col in (("NH4", 1, 2), ("NO3", 3, 4)):
            image_path = _find_representative_image(replicate_dir, condition)
            colorbar_path = replicate_dir / f"CBar_{condition}.tif"
            if not colorbar_path.exists():
                raise FileNotFoundError(colorbar_path)

            fig.add_trace(go.Image(z=tifffile.imread(image_path)), row=row, col=image_col)
            images_to_replace.append(image_path)
            _add_colorbar(fig, row, colorbar_col, colorbar_path)
            fig.update_xaxes(visible=False, showgrid=False, showline=False, row=row, col=image_col)
            fig.update_yaxes(visible=False, showgrid=False, showline=False, row=row, col=image_col)

            profile_group = df[(df["Figure"] == replicate) & (df["Treatment"] == condition)]
            if profile_group.empty:
                raise ValueError(f"No measured data for {replicate} / {condition}")
            counts[condition].append(len(profile_group))

        for condition in CONDITIONS:
            profile_group = df[(df["Figure"] == replicate) & (df["Treatment"] == condition)]
            summary = _summarize_profile(profile_group)
            x = np.arange(len(RING_COLUMNS)) * RING_DISTANCE_NM
            mean = np.asarray(summary["mean"])
            ci95 = np.asarray(summary["ci95"])
            color = CONDITION_COLORS[condition]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mean,
                    mode="lines",
                    line=dict(color=color, width=2),
                    name=condition,
                    showlegend=False,
                    hovertemplate="Distance: %{x} nm<br>Norm. intensity: %{y:.2f}<extra></extra>",
                ),
                row=row,
                col=5,
            )
            fig.add_trace(
                go.Scatter(
                    x=np.concatenate([x, x[::-1]]),
                    y=np.concatenate([mean + ci95, (mean - ci95)[::-1]]),
                    fill="toself",
                    fillcolor=("rgba(0, 106, 56, 0.18)" if condition == "NH4" else "rgba(0, 48, 135, 0.18)"),
                    line=dict(color="rgba(255,255,255,0)"),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row,
                col=5,
            )

        n_text = (
            f"<span style='color:{CONDITION_COLORS['NH4']}'>"
            f"NH<sub>4</sub> n={counts['NH4'][row - 1]}</span>, "
            f"<span style='color:{CONDITION_COLORS['NO3']}'>"
            f"NO<sub>3</sub> n={counts['NO3'][row - 1]}</span>"
        )
        fig.add_annotation(
            x=0.98,
            y=0.98,
            xref="x domain",
            yref="y domain",
            text=n_text,
            showarrow=False,
            xanchor="right",
            yanchor="top",
            borderpad=0.2,
            bgcolor="white",
            font=dict(size=PT7),
            row=row,
            col=5,
        )

        fig.update_xaxes(
            range=[0, 520],
            dtick=130,
            showgrid=False,
            row=row,
            col=5,
        )
        fig.update_yaxes(
            range=[0, 1],
            showgrid=False,
            title_text="Norm. intensity",
            title_standoff=0,
            row=row,
            col=5,
        )

    # One global legend above the panels.  Per-replicate counts are shown in
    # the corresponding profile panel.
    for condition in CONDITIONS:
        condition_label = "NH<sub>4</sub>" if condition == "NH4" else "NO<sub>3</sub>"
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines",
                line=dict(color=CONDITION_COLORS[condition], width=2),
                name=condition_label,
                legendgroup=condition,
                showlegend=True,
            ),
            row=1,
            col=5,
        )

    # Add group headers centered above each logical column.
    for xref, text in (
        ("x domain", "<b>NH<sub>4</sub></b>"),
        ("x3 domain", "<b>NO<sub>3</sub></b>"),
    ):
        fig.add_annotation(
            x=0.5,
            y=1.06,
            xref=xref,
            yref="paper",
            text=text,
            showarrow=False,
            font=dict(size=PT7 + 1),
        )

    for row, replicate in enumerate(replicates, start=1):
        yaxis = "yaxis" if row == 1 else f"yaxis{(row - 1) * 5 + 1}"
        y_domain = tuple(getattr(fig.layout, yaxis).domain)
        fig.add_annotation(
            x=-0.025,
            y=sum(y_domain) / 2,
            xref="paper",
            yref="paper",
            text=replicate,
            showarrow=False,
            xanchor="right",
            font=dict(size=PT7),
        )

    fig.update_xaxes(title_text="Distance to probe [nm]", row=len(replicates), col=5)
    fig.update_layout(
        width=680,
        height=430 if len(replicates) == 2 else 220 * len(replicates),
        margin=dict(l=68, r=8, t=76, b=48),
        legend=dict(
            orientation="h",
            x=1.0,
            xanchor="right",
            y=1.17,
            yanchor="top",
            bgcolor="rgba(0,0,0,0)",
            font=dict(size=PT7),
        ),
    )
    return fig, images_to_replace


def write_outputs(
    distribution_file: str | Path,
    output_svg: str | Path,
    output_html: str | Path,
    image_root: str | Path = f"initialData/AllFISHData/{PROBE}",
) -> None:
    output_svg = Path(output_svg)
    output_html = Path(output_html)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig, images_to_replace = build_figure(distribution_file, image_root)
    fig.write_html(output_html)
    tmp_svg = Path(f"{output_svg}.tmp")
    fig.write_image(tmp_svg, format="svg")
    replacements = replace_square_images(tmp_svg, images_to_replace, output_svg)
    tmp_svg.unlink()
    if replacements != len(images_to_replace):
        raise RuntimeError(
            f"Embedded {replacements} representative images, expected {len(images_to_replace)}"
        )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", default="Pipeline/FISH/all_distributions.tsv")
    parser.add_argument(
        "--svg",
        default="Pipeline/Paper/Supplementary/nrtBFISH/nrtBFISH.svg",
    )
    parser.add_argument(
        "--html",
        default="Pipeline/Paper/Supplementary/nrtBFISH/nrtBFISH.html",
    )
    parser.add_argument("--image-root", default=f"initialData/AllFISHData/{PROBE}")
    return parser.parse_args()


if __name__ == "__main__":
    if "snakemake" in globals():
        write_outputs(
            snakemake.input.tsv,
            snakemake.output.svg,
            snakemake.output.html,
            snakemake.params.image_root,
        )
    else:
        args = _parse_args()
        write_outputs(args.tsv, args.svg, args.html, args.image_root)
