#!/usr/bin/env python3

"""Plot representative images and per-cell 16S rRNA TAMRA intensities."""

import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import tifffile

try:
    from scripts.plotly_template import PT7, PT8, register_template
    from scripts.plotlyTrash import replace_square_images
except ModuleNotFoundError:  # Support direct execution from the scripts directory.
    from plotly_template import PT7, PT8, register_template
    from plotlyTrash import replace_square_images


register_template()

REQUIRED_COLUMNS = {
    "replicate",
    "timepoint",
    "cellidx",
    "imagename",
    "meanintensity",
}
LIGHT_CYCLE_TIMEPOINTS = ("1h_light", "6h_light", "1h_dark", "6h_dark")
TIMEPOINT_LABELS = {
    "1h_light": "1 h light",
    "6h_light": "6 h light",
    "1h_dark": "1 h dark",
    "6h_dark": "6 h dark",
    "continuous_light": "Continuous light",
}
PANEL_REPLICATES = ("Replicate 1", "Replicate 2", "Replicate 3")
PANEL_IMAGE_NUMBERS = {
    ("Replicate 1", "1h_light"): 10,
    ("Replicate 1", "6h_light"): 51,
    ("Replicate 1", "1h_dark"): 95,
    ("Replicate 1", "6h_dark"): 143,
    ("Replicate 2", "1h_light"): 18,
    ("Replicate 2", "6h_light"): 63,
    ("Replicate 2", "1h_dark"): 112,
    ("Replicate 2", "6h_dark"): 152,
    ("Replicate 3", "1h_light"): 26,
    ("Replicate 3", "6h_light"): 72,
    ("Replicate 3", "1h_dark"): 117,
    ("Replicate 3", "6h_dark"): 160,
}
PANEL_DISPLAY_MAX = {
    ("Replicate 1", "1h_light"): 4300,
    ("Replicate 1", "6h_light"): 2600,
    ("Replicate 1", "1h_dark"): 1000,
    ("Replicate 1", "6h_dark"): 1100,
    ("Replicate 2", "1h_light"): 3200,
    ("Replicate 2", "6h_light"): 4700,
    ("Replicate 2", "1h_dark"): 1300,
    ("Replicate 2", "6h_dark"): 2600,
    ("Replicate 3", "1h_light"): 7200,
    ("Replicate 3", "6h_light"): 2600,
    ("Replicate 3", "1h_dark"): 1800,
    ("Replicate 3", "6h_dark"): 1500,
}
PANEL_IMAGE_SIZE = 286
BOX_Y_MAX = 2000
ONE_POINT_PX = 4 / 3
FONT_TICK = PT8 + 3.5 + 2 * ONE_POINT_PX
FONT_PANEL_TITLE = PT8 + 3 + 2 * ONE_POINT_PX
FONT_AXIS_LABEL = PT8 + 5.5 + 2 * ONE_POINT_PX
FONT_HEADING = PT8 + 4.5 + 2 * ONE_POINT_PX
PANEL_LETTER_FONT = 18 + 2 * ONE_POINT_PX


def _replicate_number(replicate):
    match = re.search(r"(\d+)$", replicate)
    if match is None:
        raise ValueError(f"Cannot extract a replicate number from {replicate!r}")
    return int(match.group(1))


def _condition_colors(cfg):
    return {
        "1h_light": cfg["minusPColor"],
        "6h_light": cfg["membraneColor"],
        "1h_dark": cfg["cytosolColor"],
        "6h_dark": cfg["plusPColor"],
        "continuous_light": cfg["seahawksGrey"],
    }


def _validate_data(df, input_path):
    missing_columns = REQUIRED_COLUMNS.difference(df.columns)
    if missing_columns:
        raise ValueError(f"Missing columns in {input_path}: {sorted(missing_columns)}")

    unexpected_timepoints = set(df["timepoint"].unique()).difference(
        {*LIGHT_CYCLE_TIMEPOINTS, "continuous_light"}
    )
    if unexpected_timepoints:
        raise ValueError(f"Unexpected timepoints: {sorted(unexpected_timepoints)}")

    if df.empty:
        raise ValueError(f"No intensity data found in {input_path}")
    if not np.isfinite(df["meanintensity"]).all():
        raise ValueError("meanintensity contains missing or infinite values")
    if (df["meanintensity"] < 0).any():
        raise ValueError("meanintensity contains negative values")


def _collect_panel_sources(panel_paths):
    expected_by_number = {
        (_replicate_number(replicate), image_number): (replicate, timepoint)
        for (replicate, timepoint), image_number in PANEL_IMAGE_NUMBERS.items()
    }
    panel_sources = {}
    for value in panel_paths:
        path = Path(value)
        replicate_match = re.search(r"Replicate\s*(\d+)", str(path.parent))
        image_match = re.search(r"Image\s+(\d+)", path.name)
        if replicate_match is None or image_match is None:
            continue
        lookup_key = (int(replicate_match.group(1)), int(image_match.group(1)))
        panel_key = expected_by_number.get(lookup_key)
        if panel_key is None:
            continue
        if panel_key in panel_sources:
            raise ValueError(f"Duplicate representative panel for {panel_key}")
        panel_sources[panel_key] = path

    missing = sorted(set(PANEL_IMAGE_NUMBERS).difference(panel_sources))
    if missing:
        raise ValueError(f"Missing representative panels: {missing}")
    return panel_sources


def _prepare_panel_images(panel_sources, panel_dir):
    panel_dir.mkdir(parents=True, exist_ok=True)
    panels = {}
    source_rows = []
    for replicate in PANEL_REPLICATES:
        for timepoint in LIGHT_CYCLE_TIMEPOINTS:
            key = (replicate, timepoint)
            source = panel_sources[key]
            rgb = tifffile.imread(source)
            expected_shape = (PANEL_IMAGE_SIZE, PANEL_IMAGE_SIZE, 3)
            if rgb.shape != expected_shape or rgb.dtype != np.uint8:
                raise ValueError(
                    f"Representative panel {source} has shape/dtype {rgb.shape}/{rgb.dtype}; "
                    f"expected {expected_shape}/uint8"
                )

            output_name = f"replicate_{_replicate_number(replicate)}_{timepoint}.tif"
            panel_path = panel_dir / output_name
            shutil.copy2(source, panel_path)
            panels[key] = (rgb, panel_path)
            source_rows.append(
                {
                    "replicate": replicate,
                    "timepoint": timepoint,
                    "imagename": f"Image {PANEL_IMAGE_NUMBERS[key]}",
                    "display_max": PANEL_DISPLAY_MAX[key],
                    "source_image": str(source),
                    "rendered_panel": str(panel_path),
                }
            )

    pd.DataFrame(source_rows).to_csv(
        panel_dir / "PanelSources.tsv", sep="\t", index=False
    )
    return panels


def _add_boxplot(fig, values, timepoint, row, col, position, color):
    label = TIMEPOINT_LABELS[timepoint]
    fig.add_trace(
        go.Box(
            x=np.full(values.size, position),
            y=values,
            width=0.55,
            fillcolor=color,
            opacity=0.8,
            line=dict(color="black", width=0.8),
            boxpoints=False,
            whiskerwidth=0.65,
            name=label,
            showlegend=False,
            hovertemplate=(
                f"{label}<br>Mean TAMRA intensity / cell: %{{y:.1f}} A.U."
                "<extra></extra>"
            ),
        ),
        row=row,
        col=col,
    )


def _format_box_xaxis(fig, timepoints, counts, row, col):
    tickvals = list(range(len(timepoints)))
    ticktext = [
        f"{TIMEPOINT_LABELS[timepoint]}<br>n={counts[timepoint]}"
        for timepoint in timepoints
    ]
    axis_range = (
        [-2.05, 2.05]
        if len(timepoints) == 1
        else [-0.55, len(timepoints) - 0.45]
    )
    fig.update_xaxes(
        range=axis_range,
        tickmode="array",
        tickvals=tickvals,
        ticktext=ticktext,
        tickfont=dict(size=FONT_TICK),
        tickangle=-28 if len(timepoints) > 1 else 0,
        row=row,
        col=col,
    )


def _figure_specs():
    image_row = []
    for _ in LIGHT_CYCLE_TIMEPOINTS:
        image_row.extend(
            [{"type": "image", "colspan": 2}, None, {"type": "xy"}]
        )
    box_row = []
    for _ in PANEL_REPLICATES:
        box_row.extend([{"type": "xy", "colspan": 4}, None, None, None])
    blank_row = [None] * 12
    return (
        [list(image_row) for _ in PANEL_REPLICATES]
        + [blank_row, box_row, list(blank_row), list(box_row)]
    )


def _add_image_grid(fig, panels):
    image_sources = []
    image_columns = (1, 4, 7, 10)
    colorbar_columns = (3, 6, 9, 12)

    for row, replicate in enumerate(PANEL_REPLICATES, start=1):
        for timepoint, image_col, colorbar_col in zip(
            LIGHT_CYCLE_TIMEPOINTS, image_columns, colorbar_columns
        ):
            rgb, image_path = panels[(replicate, timepoint)]
            fig.add_trace(go.Image(z=rgb, hoverinfo="skip"), row=row, col=image_col)
            image_sources.append(str(image_path))
            fig.update_xaxes(visible=False, row=row, col=image_col)
            fig.update_yaxes(visible=False, row=row, col=image_col)

            display_max = PANEL_DISPLAY_MAX[(replicate, timepoint)]
            values = np.linspace(0, display_max, 256)
            fig.add_trace(
                go.Heatmap(
                    z=values.reshape(-1, 1),
                    x=[0],
                    y=values,
                    colorscale=[[0, "black"], [1, "rgb(0,255,0)"]],
                    showscale=False,
                    hovertemplate="%{y:.0f} A.U.<extra></extra>",
                ),
                row=row,
                col=colorbar_col,
            )
            colorbar_subplot = fig.get_subplot(row, colorbar_col)
            colorbar_domain = colorbar_subplot.xaxis.domain
            colorbar_width = (colorbar_domain[1] - colorbar_domain[0]) * 0.22
            fig.update_xaxes(
                domain=[colorbar_domain[0], colorbar_domain[0] + colorbar_width],
                visible=False,
                row=row,
                col=colorbar_col,
            )
            fig.update_yaxes(
                range=[0, display_max],
                tickmode="array",
                tickvals=np.linspace(0, display_max, 5),
                tickformat=".0f",
                tickfont=dict(size=FONT_TICK),
                side="right",
                ticks="outside",
                ticklen=2,
                title_text=None,
                row=row,
                col=colorbar_col,
            )

    for timepoint, image_col in zip(LIGHT_CYCLE_TIMEPOINTS, image_columns):
        fig.add_annotation(
            x=0.5,
            y=1.09,
            xref="x domain",
            yref="y domain",
            text=f"<b>{TIMEPOINT_LABELS[timepoint]}</b>",
            showarrow=False,
            xanchor="center",
            yanchor="bottom",
            font=dict(size=FONT_HEADING),
            row=1,
            col=image_col,
        )

    for row, replicate in enumerate(PANEL_REPLICATES, start=1):
        subplot = fig.get_subplot(row, 1)
        fig.add_annotation(
            x=subplot.xaxis.domain[0] - 0.017,
            y=sum(subplot.yaxis.domain) / 2,
            xref="paper",
            yref="paper",
            text=f"<i>{replicate}</i>",
            showarrow=False,
            xanchor="right",
            yanchor="middle",
            font=dict(size=FONT_HEADING),
        )

    fig.add_annotation(
        x=-0.055,
        y=1.025,
        xref="paper",
        yref="paper",
        text="<b>a</b>",
        showarrow=False,
        xanchor="left",
        yanchor="bottom",
        font=dict(size=PANEL_LETTER_FONT),
    )
    return image_sources


def _add_boxplots(fig, df, cfg):
    colors = _condition_colors(cfg)
    light_cycle = df[df["timepoint"].isin(LIGHT_CYCLE_TIMEPOINTS)]
    continuous_light = df[df["timepoint"] == "continuous_light"]
    box_columns = (1, 5, 9)

    for col, replicate in zip(box_columns, PANEL_REPLICATES):
        replicate_data = light_cycle[light_cycle["replicate"] == replicate]
        counts = {}
        for position, timepoint in enumerate(LIGHT_CYCLE_TIMEPOINTS):
            values = replicate_data.loc[
                replicate_data["timepoint"] == timepoint, "meanintensity"
            ].to_numpy(dtype=float)
            if values.size == 0:
                raise ValueError(f"No data for {replicate}, {timepoint}")
            counts[timepoint] = values.size
            _add_boxplot(
                fig,
                values,
                timepoint,
                row=5,
                col=col,
                position=position,
                color=colors[timepoint],
            )
        _format_box_xaxis(fig, LIGHT_CYCLE_TIMEPOINTS, counts, row=5, col=col)

    for col, replicate in zip(box_columns, PANEL_REPLICATES):
        values = continuous_light.loc[
            continuous_light["replicate"] == replicate, "meanintensity"
        ].to_numpy(dtype=float)
        if values.size == 0:
            raise ValueError(f"No continuous-light data for {replicate}")
        _add_boxplot(
            fig,
            values,
            "continuous_light",
            row=7,
            col=col,
            position=0,
            color=colors["continuous_light"],
        )
        _format_box_xaxis(
            fig,
            ("continuous_light",),
            {"continuous_light": values.size},
            row=7,
            col=col,
        )

    panel_letter = iter("bcdefg")
    for row, condition_label in ((5, "Light cycle"), (7, "Continuous light")):
        for col, replicate in zip(box_columns, PANEL_REPLICATES):
            fig.update_yaxes(
                range=[0, BOX_Y_MAX],
                tickmode="array",
                tickvals=[0, 500, 1000, 1500, 2000],
                tickformat="~s",
                tickfont=dict(size=FONT_TICK),
                title_text=None,
                row=row,
                col=col,
            )
            if col != box_columns[0]:
                fig.update_yaxes(showticklabels=False, row=row, col=col)
            fig.add_annotation(
                x=0.5,
                y=1.06,
                xref="x domain",
                yref="y domain",
                text=(
                    f"<b>{condition_label}</b><br>"
                    f"{replicate.lower()}"
                ),
                showarrow=False,
                xanchor="center",
                yanchor="bottom",
                font=dict(size=FONT_PANEL_TITLE),
                row=row,
                col=col,
            )
            fig.add_annotation(
                x=-0.09,
                y=1.04,
                xref="x domain",
                yref="y domain",
                text=f"<b>{next(panel_letter)}</b>",
                showarrow=False,
                xanchor="left",
                yanchor="bottom",
                font=dict(size=PANEL_LETTER_FONT),
                row=row,
                col=col,
            )

    upper_box_domain = fig.get_subplot(5, box_columns[0]).yaxis.domain
    lower_box_domain = fig.get_subplot(7, box_columns[0]).yaxis.domain
    fig.add_annotation(
        x=-0.078,
        y=(upper_box_domain[1] + lower_box_domain[0]) / 2,
        xref="paper",
        yref="paper",
        text="Mean TAMRA intensity / cell [A.U.]",
        textangle=-90,
        showarrow=False,
        xanchor="center",
        yanchor="middle",
        font=dict(size=FONT_AXIS_LABEL),
    )
    fig.add_annotation(
        x=0.5,
        y=-0.035,
        xref="paper",
        yref="paper",
        text="Lighting condition / time point",
        showarrow=False,
        xanchor="center",
        yanchor="top",
        font=dict(size=FONT_AXIS_LABEL),
    )


def build_figure(df, panels, cfg, layout):
    _validate_data(df, "intensity table")
    fig = make_subplots(
        rows=7,
        cols=12,
        specs=_figure_specs(),
        row_heights=[0.16, 0.16, 0.16, 0.04, 0.20, 0.08, 0.20],
        column_widths=[1] * 12,
        horizontal_spacing=0.005,
        vertical_spacing=0.015,
    )
    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(**layout)
    fig.update_layout(showlegend=False, hovermode="closest")

    image_sources = _add_image_grid(fig, panels)
    _add_boxplots(fig, df, cfg)
    return fig, image_sources


def main(input_tsv, panel_paths, output_svg, output_html, panel_dir, cfg, layout):
    input_tsv = Path(input_tsv)
    df = pd.read_csv(input_tsv, sep="\t")
    _validate_data(df, input_tsv)
    panel_sources = _collect_panel_sources(panel_paths)
    panels = _prepare_panel_images(panel_sources, Path(panel_dir))
    fig, images_to_replace = build_figure(df, panels, cfg, layout)

    output_svg = Path(output_svg)
    output_html = Path(output_html)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(output_html)
    temporary_svg = Path(f"{output_svg}.tmp")
    fig.write_image(
        temporary_svg,
        format="svg",
        width=fig.layout.width,
        height=fig.layout.height,
    )
    replacements = replace_square_images(temporary_svg, images_to_replace, output_svg)
    temporary_svg.unlink()
    if replacements != len(images_to_replace):
        raise RuntimeError(
            f"Expected to replace {len(images_to_replace)} raster panels in the SVG, "
            f"but replaced {replacements}."
        )


if __name__ == "__main__":
    if "snakemake" in globals():
        input_tsv = snakemake.input.tsv
        panel_paths = snakemake.input.panels
        output_svg = snakemake.output.svg
        output_html = snakemake.output.html
        panel_dir = snakemake.output.panel_images
        cfg = snakemake.config
        layout = snakemake.params.layout
    else:
        import yaml

        input_tsv = "Pipeline/FISH/16S_rRNA_light_dark_mean_intensity.tsv"
        panel_paths = sorted(
            Path("initialData/16sRNAFISHLightDark/RepresentativePanels").glob(
                "Replicate*/Composite_*.tif"
            )
        )
        with open("configFile.yml", encoding="utf-8") as handle:
            cfg = yaml.safe_load(handle)
        output_svg = cfg["Supplement16SLightDarkFISH"]["file"]
        output_html = str(Path(output_svg).with_suffix(".html"))
        panel_dir = Path(output_svg).parent / "PanelImages"
        layout = cfg["Supplement16SLightDarkFISH"]["layout"]

    main(input_tsv, panel_paths, output_svg, output_html, panel_dir, cfg, layout)
