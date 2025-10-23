import numpy as np
from scipy import stats
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import os
import tifffile
import plotly.express as px
import plotly.io as pio



def main(spots_df, cfg, outfile):
    probe = "cupA"
    treatments = ["untreated", "Puromycin"]

    # Filter dataframe for selected probe
    df = spots_df[spots_df["Probe"] == probe]
    df = df.groupby(["Treatment", "Figure"])["ProbeIntensity"].mean()

    # Create subplot grid: 2 rows, 2 columns
    fig = make_subplots(
        rows=2,
        cols=2,
        row_heights=[0.6, 0.4],  # top row for images, bottom for boxplot
        specs=[
            [{"type": "image"}, {"type": "image"}],  # Row 1: two image plots
            [{"type": "box", "colspan": 2}, None]  # Row 2: one boxplot spanning both cols
        ],
        horizontal_spacing=0.05,
        vertical_spacing=0.02,
        subplot_titles=[f"{treatments[0]}", f"{treatments[1]}", None]
    )

    # Add images for each treatment
    for i, treatment in enumerate(treatments, start=1):
        image_dir = f"initialData/AllFISHData/{probe}/Replicate 2"
        files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatment}_green.tif") and "Image" in f]
        if not files:
            print(f"No image found for {treatment} in {image_dir}")
            continue

        image_file = os.path.join(image_dir, files[0])
        fig.add_trace(
            go.Image(z=tifffile.imread(image_file)),
            row=1,
            col=i
        )
    colors = {"untreated": cfg["membraneColor"], "Puromycin": cfg["cytosolColor"]}
    # Add boxplot spanning both columns
    for treatment in treatments:
        fig.add_trace(
            go.Box(
                y=df[treatment],
                name=treatment,
                boxmean=True,
                marker=dict(color=colors[treatment]),
            ),
            row=2,
            col=1
        )
    fig.update_xaxes(showticklabels=False, row=1)
    fig.update_yaxes(showticklabels=False, row=1)
    fig.update_yaxes(title=dict(text="mean Intensities"), row=2)
    fig.update_layout(cfg["LayoutDefault"])

    fig.update_layout(cfg["F4CupAIntensities"]["layout"])
    fig.update_annotations(cfg["F4CupAIntensities"]["annotations"])
    fig.write_image(outfile)



if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        cfg = snakemake.config
        outfile = snakemake.output.svg
    else:
        import yaml

        file = "Pipeline/FISH/single_spots.tsv"
        outfile = "foo.svg"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)


    spots_df = pd.read_csv(file, sep="\t")
    main(spots_df, cfg, outfile)
