import numpy as np
from scipy import stats
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import os
import tifffile
import plotly.express as px
import plotly.io as pio
from scipy.stats import ttest_rel, ttest_ind
from statsmodels.stats.multitest import multipletests
from scripts.plotly_template import register_template, PT7
import pytesseract

register_template()

pio.renderers.default = "browser"

def _color_to_calpha(color: str, alpha: float = 0.2):
    color = color.split("(")[-1].split(")")[0]
    return f"rgba({color}, {alpha})"


def find_and_extract_cbar(image_dir, treatment):
    cbar_files = [f for f in os.listdir(image_dir) if "CBar" in f]
    if len(cbar_files) > 1:
        cbarfile = [f for f in cbar_files if treatment in f][0]
    else:
        cbarfile = cbar_files[0]
    image_file = os.path.join(image_dir, cbarfile)
    img = tifffile.imread(image_file)
    max_v = int(pytesseract.image_to_string(img).split("\n")[0])
    mask = (img[..., 0] != 255) & (img[..., 2] != 255)
    max_color = img[mask, 1].max()
    fake_heat = np.linspace(0, max_v, 100).reshape(100, 1)
    trace = go.Heatmap(
                z=fake_heat,
                colorscale=["rgb(0, 0, 0)", f"rgb(0, {max_color}, 0)"],
                showscale=False,
                # Hide the actual heatmap by making it tiny and offset
                y=fake_heat.flatten(),  # off the visible range
                x=[-1],
                # Configure colorbar
            )
    return trace, max_v

def main(results_df, cfg, outfile, outtsv):
    probes = ["psbA2", "cupA"]
    treatments = ["Untreated", "Puromycin"]
    # Create subplot grid: 2 rows, 2 columns
    r = 0.08
    l = 0.001
    base = 1/len(probes)
    column_widths = [0.4 * base, 0.4 * base, 0.2 *base] * len(probes)

    fig = make_subplots(
        rows=2,
        cols=6,
        column_widths=column_widths,
        row_heights=[0.6, 0.4],  # top row for images, bottom for boxplot
        specs=[
            [{"type": "image", "l": l, "r": l}, {"type": "image", "l": l, "r": l}, {"type": "xy", "b": 0.02, "t": 0.05, "r": r}] * len(probes) ,
            [{"type": "box", "colspan": 3, "r":(r+l )}, None,  None] * 2  # Row 2: one boxplot spanning both cols
        ],
        horizontal_spacing=0,
        vertical_spacing=0.05,
        subplot_titles=[f"Untreated", "Puromycin", None] * 2
    )
    fig.update_annotations(font=dict(size=PT7))

    # Add images for each treatment
    idx = 1
    for probe in probes:
        for i, treatment in enumerate(treatments, start=1):
            treatl = "untreated" if treatment == "Untreated" else treatment
            image_dir = f"initialData/AllFISHData/{probe}/Replicate 3"
            files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatl}_green.tif") and "Image" in f]
            if not files:
                print(f"No image found for {treatment} in {image_dir}")
                continue


            image_file = os.path.join(image_dir, files[0])
            fig.add_trace(
                go.Image(z=tifffile.imread(image_file)),
                row=1,
                col=idx
            )
            fig.update_xaxes(
                ticks="",
                showticklabels=False,
                row=1,
                col=idx
            )
            fig.update_yaxes(
                ticks="",
                showticklabels=False,
                row=1,
                col=idx,
                showline=False
            )
            if treatment == "Puromycin":
                colorscale, max_v = find_and_extract_cbar(image_dir, treatment)
                fig.add_trace(colorscale, row=1, col=idx+1)
                fig.update_yaxes(
                    tickvals=[0, max_v // 2, max_v], row=1, col=idx+1,
                    range=[0, max_v],
                    side="right",
                    showticklabels=True,
                    showline=False

                )
                idx += 1
            idx += 1

    colors = {"Untreated": cfg["minusPColor"], "Puromycin": cfg["plusPColor"]}
    # Add boxplot spanning both columns
    for pidx , probe in enumerate(probes):
        for treatment in treatments:
            intensities = eval(results_df[results_df.Probe == probe][treatment].iloc[0])
            fig.add_trace(
                go.Box(
                    y=intensities,
                    name=f"{treatment} – <i>{probe}</i>",
                    boxmean=True,
                    fillcolor=_color_to_calpha(colors[treatment], 0.7),
                    line=dict(width=1, color="black")
                ),
                row=2,
                col=1 if pidx == 0 else 4
            )
            if treatment == "Puromycin":
                pval = results_df.loc[results_df["Probe"] == probe, "adj p-value (FDR)"].values
                if len(pval) > 0 and pval[0] < 0.05:
                    # Get approximate y position (top of box)
                    y_max = max(intensities) * 1.2
                    # choose symbol by significance level
                    p = pval[0]
                    if p < 0.0001:
                        stars = "****"
                    elif p < 0.001:
                        stars = "***"
                    elif p < 0.01:
                        stars = "**"
                    elif p < 0.05:
                        stars = "*"
                    else:
                        stars = "n.s."

                    # Add a star annotation above the box
                    fig.add_annotation(
                        text=stars,  # or "*" if you prefer plain text
                        x=f"{treatment} – <i>{probe}</i>",
                        y=y_max,
                        xref="x",
                        yref="y",
                        showarrow=False,
                        font=dict(size=10, color="black"),
                        row=2,
                        col=1 if pidx == 0 else 4
                    )
    fig.update_yaxes(title=dict(text="Mean intensities"), row=2, col=1)


    hs = 0.1
    y = 1.1

    fig.add_shape(
        type="line",
        xref="x7",
        yref="paper",
        x0=0,
        x1=1,
        y0=y,
        y1=y,
        line=dict(color="black", width=1),
        opacity=1
    )
    fig.add_shape(
        type="line",
        xref="x8",
        yref="paper",
        x0=0,
        x1=1,
        y0=y,
        y1=y,
        line=dict(color="black", width=1),
        opacity=1
    )
    fig.add_annotation(
        xref="x7",
        yref="paper",
        x=0.5,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=f"<i>{probes[0]}</i>",
        showarrow=False,
    )
    fig.add_annotation(
        xref="x8",
        yref="paper",
        x=0.5,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=f"<i>{probes[1]}</i>",
        showarrow=False,
    )



    fig.update_layout(cfg["LayoutDefault"])

    fig.update_xaxes(ticks="", row=1, showticklabels=False, showline=False)
    fig.update_layout(cfg["F4CupAIntensities"]["layout"])
    fig.write_image(outfile)
    results_df.to_csv(outtsv, index=False, sep="\t")




if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        cfg = snakemake.config
        outfile = snakemake.output.svg
        outtsv = snakemake.output.tsv
    else:
        import yaml

        file = "Pipeline/Paper/Supplementary/FISHIntensitiesSummary.tsv"
        outfile = "foo.svg"
        outtsv = "foo.tsv"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)


    spots_df = pd.read_csv(file, sep="\t")
    main(spots_df, cfg, outfile, outtsv)
