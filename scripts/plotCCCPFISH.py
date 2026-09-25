#!/usr/bin/env python

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import tifffile
from skimage.measure import regionprops_table

try:
    from scripts.plotly_template import PT7, PT8, register_template
    from scripts.plotlyTrash import replace_images
except ModuleNotFoundError:
    from plotly_template import PT7, PT8, register_template
    from plotlyTrash import replace_images


PROBES = ("psaA", "psbA2", "cupA")
PROBE_LABELS = {"psaA": "psaA", "psbA2": "psbA", "cupA": "cupA"}
TREATMENTS = ("DMSO", "CCCP")
REPLICATE = "Replicate 4"
DISTANCE_STEP_NM = 26
IMAGE_VERTICAL_SPACING = 0.025
COLORBAR_VERTICAL_INSET = 0.015
PANEL_IMAGE_FILES = {
    ("psaA", "DMSO"): "Composite_Image 28_Out-Channel Alignment-10_0.553-1_DMSO_psaA.tif",
    ("psaA", "CCCP"): "Composite_Image 9_Out-Channel Alignment-02_0.557-1_CCCP_psaA.tif",
    ("psbA2", "DMSO"): "Composite_Image 20_Out-Channel Alignment-01_0.561-1_DMSO_psbA2.tif",
    ("psbA2", "CCCP"): "Composite_Image 57_Out-Channel Alignment-10_0.52-1_CCCP_psbA2.tif",
    ("cupA", "DMSO"): "Composite_Image 13_Out-Channel Alignment-07_0.561-1_DMSO_cupA.tif",
    ("cupA", "CCCP"): "Composite_Image 67_Out-Channel Alignment-07_0.52-1_CCCP_cupA.tif",
}
COLORBAR_FILES = {probe: f"CBar_CCCP_{probe}.tif" for probe in PROBES}

DMSO_COLOR = "rgb(0, 106, 56)"
CCCP_COLOR = "rgb(0, 48, 135)"
VIOLIN_COLORS = {
    "DMSO": "rgba(0, 106, 56, 0.82)",
    "CCCP": "rgba(0, 48, 135, 0.96)",
}


def _median_intensity(mask, intensity_image):
    return np.median(intensity_image[mask])


def _image_key(path):
    return path.name.split("_")[1]


def _collect_image_records(image_paths):
    records = {}
    for value in image_paths:
        path = Path(value)
        if path.suffix.lower() not in {".tif", ".tiff"}:
            continue
        parts = path.parts
        try:
            fish_index = parts.index("AllFISHData")
        except ValueError:
            continue

        probe = parts[fish_index + 1]
        if probe not in PROBES or REPLICATE not in parts:
            continue
        treatment_dir = path.parent.name
        if not treatment_dir.startswith("Segmentation_"):
            continue
        treatment = treatment_dir.removeprefix("Segmentation_")
        if treatment not in TREATMENTS:
            continue

        name = path.name
        if "_corrected_total_cell" in name:
            kind = "TotalCell"
        elif "_label_image_corrected" in name:
            kind = "Probe"
        elif name.startswith("Composite_Image"):
            kind = "Intensity"
        else:
            continue

        key = _image_key(path)
        records.setdefault((probe, treatment, key), {})[kind] = path

    complete = {
        key: value
        for key, value in records.items()
        if {"Intensity", "Probe", "TotalCell"}.issubset(value)
    }
    expected_groups = {(probe, treatment) for probe in PROBES for treatment in TREATMENTS}
    found_groups = {(probe, treatment) for probe, treatment, _ in complete}
    missing = sorted(expected_groups - found_groups)
    if missing:
        raise ValueError(f"No complete image triplets found for: {missing}")
    return complete


def _extract_violin_data(records):
    rows = []
    for (probe, treatment, image_key), files in sorted(records.items()):
        spot_labels = tifffile.imread(files["Probe"])
        cell_labels = tifffile.imread(files["TotalCell"])
        intensity_stack = tifffile.imread(files["Intensity"])
        tamra = intensity_stack[1] if intensity_stack.ndim > 2 else intensity_stack

        cell_properties = regionprops_table(
            spot_labels,
            intensity_image=cell_labels,
            properties=("label",),
            extra_properties=(_median_intensity,),
        )
        intensity_properties = regionprops_table(
            spot_labels,
            intensity_image=tamra,
            properties=("label",),
            extra_properties=(_median_intensity,),
        )
        spots = pd.DataFrame(cell_properties).rename(
            columns={"_median_intensity": "CellIdx"}
        )
        spot_intensities = pd.DataFrame(intensity_properties).rename(
            columns={"_median_intensity": "ProbeIntensity"}
        )
        spots = spots.merge(spot_intensities, on="label", how="inner")
        spots["CellIdx"] = spots["CellIdx"].astype(int)
        spots = spots[spots["CellIdx"] != 0].copy()
        background = float(np.median(tamra[spot_labels == 0]))
        spots["Value"] = spots["ProbeIntensity"] - background

        counts = spots["CellIdx"].value_counts()
        for cell_idx in np.delete(np.unique(cell_labels), 0):
            rows.append(
                {
                    "Metric": "Spots per cell",
                    "Probe": probe,
                    "Treatment": treatment,
                    "Replicate": REPLICATE,
                    "Image": image_key,
                    "CellIdx": int(cell_idx),
                    "SpotLabel": np.nan,
                    "Value": int(counts.get(int(cell_idx), 0)),
                    "SourceImage": str(files["Intensity"]),
                }
            )

        for spot in spots.itertuples(index=False):
            rows.append(
                {
                    "Metric": "TAMRA intensity per spot",
                    "Probe": probe,
                    "Treatment": treatment,
                    "Replicate": REPLICATE,
                    "Image": image_key,
                    "CellIdx": int(spot.CellIdx),
                    "SpotLabel": int(spot.label),
                    "Value": float(spot.Value),
                    "SourceImage": str(files["Intensity"]),
                }
            )
    return pd.DataFrame(rows)


def _prepare_panel_images(representative_paths, panel_dir):
    panel_dir.mkdir(parents=True, exist_ok=True)
    paths_by_name = {Path(path).name: Path(path) for path in representative_paths}
    panels = {}
    colorbars = {}
    sources = []
    for probe in PROBES:
        for treatment in TREATMENTS:
            filename = PANEL_IMAGE_FILES[(probe, treatment)]
            if filename not in paths_by_name:
                raise FileNotFoundError(f"Representative panel image is missing: {filename}")
            source = paths_by_name[filename]
            panels[(probe, treatment)] = (tifffile.imread(source), source)
            sources.append(
                {
                    "Probe": probe,
                    "Treatment": treatment,
                    "SourceImage": str(source),
                }
            )
        colorbar_filename = COLORBAR_FILES[probe]
        if colorbar_filename not in paths_by_name:
            raise FileNotFoundError(f"Representative colorbar is missing: {colorbar_filename}")
        colorbar_path = paths_by_name[colorbar_filename]
        colorbars[probe] = (tifffile.imread(colorbar_path), colorbar_path)
    pd.DataFrame(sources).to_csv(panel_dir / "PanelSources.tsv", sep="\t", index=False)
    return panels, colorbars


def _distance_columns(distribution):
    return sorted(
        (column for column in distribution.columns if column.endswith(" nm")),
        key=lambda column: int(column.split()[0]),
    )


def _rounded_upper(values, step, minimum):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return minimum
    robust_max = np.quantile(finite, 0.995)
    return max(minimum, float(np.ceil(robust_max / step) * step))


def _stacked_panel_domains(count, spacing):
    panel_height = (1 - spacing * (count - 1)) / count
    return [
        (
            1 - (index + 1) * panel_height - index * spacing,
            1 - index * (panel_height + spacing),
        )
        for index in range(count)
    ]


def _add_violin_panel(fig, data, metric, row, col, y_title, rng, upper):
    positions = {
        (probe, treatment): probe_idx * 2 + treatment_idx
        for probe_idx, probe in enumerate(PROBES)
        for treatment_idx, treatment in enumerate(TREATMENTS)
    }
    for probe in PROBES:
        probe_label = PROBE_LABELS[probe]
        for treatment in TREATMENTS:
            position = positions[(probe, treatment)]
            values = data.loc[
                (data["Metric"] == metric)
                & (data["Probe"] == probe)
                & (data["Treatment"] == treatment),
                "Value",
            ].dropna().to_numpy(dtype=float)
            fig.add_trace(
                go.Violin(
                    x=np.full(values.size, position),
                    y=values,
                    width=0.72,
                    spanmode="hard",
                    scalemode="width",
                    fillcolor=VIOLIN_COLORS[treatment],
                    line=dict(color="black", width=0.8),
                    box_visible=True,
                    meanline_visible=False,
                    points=False,
                    hovertemplate=(
                        f"<i>{probe_label}</i><br>{treatment}<br>{y_title}: %{{y:.1f}}<extra></extra>"
                    ),
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

            max_points = 700
            if values.size > max_points:
                values_for_points = rng.choice(values, max_points, replace=False)
            else:
                values_for_points = values
            jitter = rng.uniform(-0.12, 0.12, values_for_points.size)
            fig.add_trace(
                go.Scatter(
                    x=position + jitter,
                    y=values_for_points,
                    mode="markers",
                    marker=dict(color="black", size=1.5, opacity=0.55),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

    fig.update_xaxes(
        range=[-0.5, 5.5],
        tickmode="array",
        tickvals=list(range(6)),
        ticktext=["DMSO", "CCCP"] * 3,
        tickangle=-45,
        tickfont=dict(size=PT7),
        row=row,
        col=col,
    )
    fig.update_yaxes(
        range=[0, upper * 1.12],
        title_text=y_title,
        title_font=dict(size=PT7),
        title_standoff=3,
        tickfont=dict(size=PT7),
        row=row,
        col=col,
    )
    for separator in (1.5, 3.5):
        fig.add_vline(
            x=separator,
            line=dict(color="rgb(220, 220, 220)", width=0.7),
            row=row,
            col=col,
        )
    for probe_idx, probe in enumerate(PROBES):
        fig.add_annotation(
            x=probe_idx * 2 + 0.5,
            y=upper * 1.07,
            text=f"<b><i>{PROBE_LABELS[probe]}</i></b>",
            showarrow=False,
            xanchor="center",
            yanchor="middle",
            row=row,
            col=col,
        )


def build_figure(distribution, violin_data, panels, colorbars, cfg, layout):
    specs = [
        [
            {"type": "image", "rowspan": 2},
            {"type": "image", "rowspan": 2},
            {"type": "image", "rowspan": 2},
            None,
            {"type": "image", "rowspan": 2},
            None,
            {"type": "xy", "rowspan": 3, "colspan": 2},
            None,
        ],
        [None, None, None, None, None, None, None, None],
        [
            {"type": "image", "rowspan": 2},
            {"type": "image", "rowspan": 2},
            {"type": "image", "rowspan": 2},
            None,
            {"type": "xy", "rowspan": 2},
            None,
            None,
            None,
        ],
        [None, None, None, None, None, None,
         {"type": "xy", "rowspan": 3, "colspan": 2}, None],
        [
            {"type": "image", "rowspan": 2},
            {"type": "image", "rowspan": 2},
            {"type": "xy", "rowspan": 2},
            None,
            {"type": "xy", "rowspan": 2},
            None,
            None,
            None,
        ],
        [None, None, None, None, None, None, None, None],
    ]
    fig = make_subplots(
        rows=6,
        cols=8,
        specs=specs,
        row_heights=[1 / 6] * 6,
        column_widths=[0.145, 0.145, 0.055, 0.07, 0.21, 0.055, 0.14, 0.14],
        horizontal_spacing=0.02,
        vertical_spacing=0.15,
        column_titles=["DMSO", "CCCP", None, None, None, None, None, None],
    )
    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(**layout)
    fig.update_annotations(font=dict(size=PT8))

    row_starts = (1, 3, 5)
    image_domains = _stacked_panel_domains(
        len(row_starts), IMAGE_VERTICAL_SPACING
    )
    for row, domain in zip(row_starts, image_domains):
        for col in (1, 2):
            fig.update_yaxes(domain=domain, row=row, col=col)
        colorbar_domain = (
            domain[0] + COLORBAR_VERTICAL_INSET,
            domain[1] - COLORBAR_VERTICAL_INSET,
        )
        fig.update_yaxes(domain=colorbar_domain, row=row, col=3)

    distance_columns = _distance_columns(distribution)
    distance_x = np.array([int(column.split()[0]) for column in distance_columns])
    image_sources = []

    for probe, row in zip(PROBES, row_starts):
        probe_label = PROBE_LABELS[probe]
        for treatment, col in zip(TREATMENTS, (1, 2)):
            rgb, image_path = panels[(probe, treatment)]
            fig.add_trace(go.Image(z=rgb, hoverinfo="skip"), row=row, col=col)
            image_sources.append(str(image_path))
            fig.update_xaxes(visible=False, row=row, col=col)
            fig.update_yaxes(visible=False, row=row, col=col)

        colorbar, colorbar_path = colorbars[probe]
        fig.add_trace(
            go.Image(z=colorbar, hoverinfo="skip"),
            row=row,
            col=3,
        )
        image_sources.append(str(colorbar_path))
        fig.update_xaxes(visible=False, row=row, col=3)
        fig.update_yaxes(visible=False, row=row, col=3)

        probe_distribution = distribution[
            (distribution["Probe"] == probe)
            & (distribution["Replicate"] == REPLICATE)
            & (distribution["Treatment"].isin(TREATMENTS))
        ]
        peaks = {}
        sample_sizes = {}
        for treatment in TREATMENTS:
            current = probe_distribution[probe_distribution["Treatment"] == treatment]
            if current.empty:
                raise ValueError(f"Missing distance summary for {probe}, {treatment}")
            current = current.iloc[0]
            color = DMSO_COLOR if treatment == "DMSO" else CCCP_COLOR
            dash = "solid" if treatment == "DMSO" else "dash"
            fig.add_trace(
                go.Scatter(
                    x=distance_x,
                    y=current[distance_columns].to_numpy(dtype=float),
                    mode="lines",
                    line=dict(color=color, width=1.8, dash=dash),
                    hovertemplate=(
                        f"<i>{probe_label}</i><br>{treatment}<br>Distance: %{{x}} nm"
                        "<br>Normalized intensity: %{y:.3f}<extra></extra>"
                    ),
                    name=treatment,
                    legendgroup=treatment,
                    showlegend=probe == PROBES[0],
                ),
                row=row,
                col=5,
            )
            peaks[treatment] = int(current["Peak position"]) * DISTANCE_STEP_NM
            sample_sizes[treatment] = int(current["n_spots"])

        n_text = "<b>n: </b>" + ", ".join(
            f"<span style='color:{DMSO_COLOR if treatment == 'DMSO' else CCCP_COLOR}'>"
            f"{sample_sizes[treatment]}</span>"
            for treatment in TREATMENTS
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

        for peak_idx, peak in enumerate(dict.fromkeys(peaks.values())):
            fig.add_vline(
                x=peak,
                line=dict(color="black", width=0.8, dash="dot"),
                row=row,
                col=5,
            )
            fig.add_annotation(
                x=peak,
                y=0.035 + peak_idx * 0.07,
                text=str(peak),
                showarrow=False,
                xanchor="left",
                yanchor="bottom",
                font=dict(size=PT7),
                row=row,
                col=5,
            )

        fig.update_xaxes(
            range=[-20, 520],
            tickmode="array",
            tickvals=[0, 250, 500],
            title_standoff=2,
            title_text="Distance to probe [nm]",
            title_font=dict(size=PT7),
            tickfont=dict(size=PT7),
            row=row,
            col=5,
        )
        fig.update_yaxes(
            range=[0, 1],
            tickmode="array",
            tickvals=[0, 0.25, 0.5, 0.75, 1],
            title_text=f"Norm. intensity",
            tickfont=dict(size=PT7),
            row=row,
            col=5,
        )
        fig.add_annotation(
            x=0.5,
            y=0.98,
            xref="x domain",
            yref="y domain",
            text=f"<b><i>{probe_label}</i></b>",
            showarrow=False,
            xanchor="center",
            yanchor="bottom",
            row=row,
            col=5,
        )

    for row, probe in zip(row_starts, PROBES):
        image_subplot = fig.get_subplot(row, 1)
        y_domain = image_subplot.yaxis.domain
        y = sum(y_domain) / 2
        fig.add_annotation(
            x=-0.025,
            y=y,
            xref="paper",
        yref="paper",
        text=f"<i>{PROBE_LABELS[probe]}</i>",
        showarrow=False,
        textangle=-90,
        xanchor="center",
        yanchor="middle",
        font=dict(size=PT8 + 1),
        )

    spots_per_cell = violin_data.loc[
        violin_data["Metric"] == "Spots per cell", "Value"
    ]
    spot_intensities = violin_data.loc[
        violin_data["Metric"] == "TAMRA intensity per spot", "Value"
    ]
    spots_upper = _rounded_upper(spots_per_cell, step=5, minimum=40)
    intensity_upper = _rounded_upper(spot_intensities, step=500, minimum=1000)
    rng = np.random.default_rng(20260818)
    _add_violin_panel(
        fig,
        violin_data,
        "Spots per cell",
        row=1,
        col=7,
        y_title="Spots per cell",
        rng=rng,
        upper=spots_upper,
    )
    _add_violin_panel(
        fig,
        violin_data,
        "TAMRA intensity per spot",
        row=4,
        col=7,
        y_title="TAMRA intensity / spot [A.U.]",
        rng=rng,
        upper=intensity_upper,
    )
    fig.update_yaxes(domain=[0.58, 0.91], row=1, col=7)
    fig.update_yaxes(domain=[0.07, 0.40], row=4, col=7)
    line_x_domain = fig.get_subplot(1, 5).xaxis.domain
    fig.update_layout(
        showlegend=True,
        hovermode="closest",
        violinmode="overlay",
        legend=dict(
            visible=True,
            orientation="h",
            x=sum(line_x_domain) / 2,
            xanchor="center",
            y=1.015,
            yanchor="bottom",
            font=dict(size=PT7),
            bgcolor="rgba(255,255,255,0)",
            traceorder="normal",
        ),
    )
    
    fig.add_annotation(
                x=-0.02,
                y=1.07,
                xref="paper",
                yref="paper",
                text=f"<b>a</b>",
                showarrow=False,
                xanchor="center",
                yanchor="bottom",
                font=dict(size=14),

            )
    fig.add_annotation(
                    x=0.36,
                    y=1.07,
                    xref="paper",
                    yref="paper",
                    text=f"<b>b</b>",
                    showarrow=False,
                    xanchor="center",
                    yanchor="bottom",
                    font=dict(size=14),
    
                )
    fig.add_annotation(
                    x=0.66,
                    y=0.93,
                    xref="paper",
                    yref="paper",
                    text=f"<b>c</b>",
                    showarrow=False,
                    xanchor="center",
                    yanchor="bottom",
                    font=dict(size=14),
    
                )
    fig.add_annotation(
                        x=0.66,
                        y=0.44,
                        xref="paper",
                        yref="paper",
                        text=f"<b>d</b>",
                        showarrow=False,
                        xanchor="center",
                        yanchor="bottom",
                        font=dict(size=14),
        
                    )
    return fig, image_sources


def main(
    distribution_file,
    image_paths,
    representative_paths,
    outfile,
    html_outfile,
    plot_data_out,
    panel_dir,
    cfg,
    layout,
):
    register_template()
    distribution = pd.read_csv(distribution_file, sep="\t")
    records = _collect_image_records(image_paths)
    violin_data = _extract_violin_data(records)
    plot_data_out = Path(plot_data_out)
    plot_data_out.parent.mkdir(parents=True, exist_ok=True)
    violin_data.to_csv(plot_data_out, sep="\t", index=False)

    panels, colorbars = _prepare_panel_images(representative_paths, Path(panel_dir))
    fig, images_to_replace = build_figure(
        distribution, violin_data, panels, colorbars, cfg, layout
    )

    outfile = Path(outfile)
    html_outfile = Path(html_outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    html_outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(html_outfile)
    tmp_out = Path(f"{outfile}.tmp")
    fig.write_image(tmp_out, format="svg")
    replacements = replace_images(tmp_out, images_to_replace, outfile)
    tmp_out.unlink()
    if replacements != len(images_to_replace):
        raise RuntimeError(
            f"Expected to replace {len(images_to_replace)} raster panels in the SVG, "
            f"but replaced {replacements}."
        )


if __name__ == "__main__":
    if "snakemake" in globals():
        distribution_file = snakemake.input.distribution
        image_paths = snakemake.input.images
        representative_paths = snakemake.input.representatives
        outfile = snakemake.output.svg
        html_outfile = snakemake.output.html
        plot_data_out = snakemake.output.plot_data
        panel_dir = snakemake.output.panel_images
        cfg = snakemake.config
        layout = snakemake.params.layout
    else:
        import yaml

        distribution_file = "Pipeline/Paper/Supplementary/FISHDistributionSummary.tsv"
        image_paths = [
            path
            for path in Path("initialData/AllFISHData").glob(
                "*/Replicate 4/Segmentation_*/*"
            )
            if path.suffix.lower() in {".tif", ".tiff"}
        ]
        representative_paths = [
            Path("initialData/AllFISHData") / probe / REPLICATE / filename
            for (probe, _), filename in PANEL_IMAGE_FILES.items()
        ] + [
            Path("initialData/AllFISHData") / probe / REPLICATE / filename
            for probe, filename in COLORBAR_FILES.items()
        ]
        outfile = "Pipeline/Paper/Supplementary/CCCPFISH/CCCPFISH.svg"
        html_outfile = "Pipeline/Paper/Supplementary/CCCPFISH/CCCPFISH.html"
        plot_data_out = "Pipeline/Paper/Supplementary/CCCPFISH/PlotData.tsv"
        panel_dir = "Pipeline/Paper/Supplementary/CCCPFISH/PanelImages"
        with open("configFile.yml", encoding="utf-8") as handle:
            cfg = yaml.safe_load(handle)
        layout = cfg["SupplementCCCPFISH"]["layout"]

    main(
        distribution_file,
        image_paths,
        representative_paths,
        outfile,
        html_outfile,
        plot_data_out,
        panel_dir,
        cfg,
        layout,
    )
