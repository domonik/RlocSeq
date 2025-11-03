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


pio.renderers.default = "browser"


def calc_summary_and_pval(spots_df, probes):
    spots_df = spots_df[spots_df.Treatment != "Kasugamycin"]
    spots_df["subtractedIntensity"] = spots_df["ProbeIntensity"] - spots_df["background"]

    # Filter dataframe for selected probe
    df = spots_df[spots_df["Probe"].isin(probes)]
    df = df.groupby(["Treatment", "Probe", "Figure", ])["subtractedIntensity"].mean().reset_index()
    # Run Welch’s t-test for each probe
    results = []
    for probe in df['Probe'].unique():
        untreated = df.query("Treatment == 'untreated' and Probe == @probe")['subtractedIntensity']
        puromycin = df.query("Treatment == 'Puromycin' and Probe == @probe")['subtractedIntensity']
        if len(untreated) == len(puromycin):
            tstat, pval = ttest_rel(untreated, puromycin, )
        else:
            tstat, pval = ttest_ind(untreated, puromycin, equal_var=False)
        results.append({'Probe': probe, 't-statistic': tstat, 'p-value': pval, "untreated": untreated, "puromycin": puromycin})

    results_df = pd.DataFrame(results)
    reject, pvals_corr, _, _ = multipletests(results_df['p-value'], method='fdr_bh')

    results_df['adj p-value (FDR)'] = pvals_corr
    results_df['significant (FDR<0.05)'] = reject
    return df, results_df




def main(spots_df, cfg, outfile, outtsv):
    probes = ["psbA2", "cupA"]
    treatments = ["untreated", "Puromycin"]
    df, results_df = calc_summary_and_pval(spots_df, probes)
    # Create subplot grid: 2 rows, 2 columns
    fig = make_subplots(
        rows=2,
        cols=4,
        row_heights=[0.6, 0.4],  # top row for images, bottom for boxplot
        specs=[
            [{"type": "image"}, {"type": "image"}, {"type": "image"}, {"type": "image"}],  # Row 1: two image plots
            [{"type": "box", "colspan": 2}, None, {"type": "box", "colspan": 2}, None]  # Row 2: one boxplot spanning both cols
        ],
        horizontal_spacing=0.05,
        vertical_spacing=0.02,
        subplot_titles=[f"untreated", "Puromycin", "untreated", "Puromycin"]
    )

    # Add images for each treatment
    idx = 0
    for probe in probes:
        for i, treatment in enumerate(treatments, start=1):
            idx += 1
            image_dir = f"initialData/AllFISHData/{probe}/Replicate 3"
            files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatment}_green.tif") and "Image" in f]
            if not files:
                print(f"No image found for {treatment} in {image_dir}")
                continue

            image_file = os.path.join(image_dir, files[0])
            fig.add_trace(
                go.Image(z=tifffile.imread(image_file)),
                row=1,
                col=idx
            )
    colors = {"untreated": cfg["cytosolColor"], "Puromycin": cfg["membraneColor"]}
    # Add boxplot spanning both columns
    for pidx , probe in enumerate(probes):
        for treatment in treatments:
            group = df[(df.Treatment == treatment) & (df.Probe == probe)]
            name = "untreated" if treatment == "untreated" else "Puromycin"
            name = name + f" {probe}"
            fig.add_trace(
                go.Box(
                    y=group["subtractedIntensity"].tolist(),
                    name=f"{treatment} – <i>{probe}</i>",
                    boxmean=True,
                    marker=dict(color=colors.get(treatment, "black")),
                ),
                row=2,
                col=1 if pidx == 0 else 3
            )
    fig.update_xaxes(showticklabels=False, row=1)
    fig.update_yaxes(showticklabels=False, row=1)
    fig.update_yaxes(title=dict(text="mean Intensities"), row=2, col=1)


    hs = 0.1
    l = (1- (len(probes)-1)*hs) / len(probes)
    y = 1.1

    fig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=0 + hs,
        x1=0.5 - hs,
        y0=y,
        y1=y,
        line=dict(color="black", width=2),
        opacity=1
    )
    fig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=0.5 + hs,
        x1=1- hs,
        y0=y,
        y1=y,
        line=dict(color="black", width=2),
        opacity=1
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=1 / 4,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=f"<i>{probes[0]}</i>",
        showarrow=False,
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=3 / 4,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=f"<i>{probes[1]}</i>",
        showarrow=False,
    )



    fig.update_layout(cfg["LayoutDefault"])


    fig.update_layout(cfg["F4CupAIntensities"]["layout"])
    fig.update_annotations(cfg["F4CupAIntensities"]["annotations"])
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

        file = "Pipeline/FISH/single_spots.tsv"
        outfile = "foo.svg"
        outtsv = "foo.tsv"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)


    spots_df = pd.read_csv(file, sep="\t")
    main(spots_df, cfg, outfile, outtsv)
