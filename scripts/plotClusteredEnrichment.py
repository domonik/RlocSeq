from DEplots.enrichment import enrichment_plot_from_cp_table
import pandas as pd
import plotly.graph_objects as go
from scripts.plotly_template import register_template


register_template()

def main(tsv, outfile, layout, config):
    df = pd.read_csv(tsv, sep="\t")
    df["Description"] = df["Description"].str.replace("bounded", "bounded<br>", regex=False)
    fig: go.Figure = enrichment_plot_from_cp_table(df, colorscale=[config["plusPColor"], config["minusPColor"],])
    term = "outer membrane-bounded periplasmic space"
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(layout)
    fig.write_image(outfile)


if __name__ == '__main__':
    main(
        tsv=snakemake.input.tsv,
        outfile=snakemake.output.svg,
        layout=snakemake.params.layout,
        config=snakemake.config
    )