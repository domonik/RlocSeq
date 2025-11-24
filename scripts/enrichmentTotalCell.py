import pandas as pd
import numpy as np
from DEplots.enrichment import enrichment_plot_from_cp_table
from scripts.plotly_template import register_template

register_template()

def main(tsv, config, svg, html):
    colorscale = ["purple", "yellow"]

    df = pd.read_csv(tsv, sep="\t", )
    df["p.adjust"] = -1 * np.log10(df["p.adjust"])

    fig = enrichment_plot_from_cp_table(
        df,
        mode="scatter",
        colorscale=colorscale
    )

    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config["F4ClusteredGO"]["layout"])
    fig.update_layout(coloraxis=dict(colorbar=dict(title=dict(text="-log<sub>10</sub>(p.adjust)"))))


    fig.write_image(svg)
    fig.write_html(html)
    p = 0

if __name__ == '__main__':
    main(snakemake.input.tsv, config=snakemake.config, svg=snakemake.output.svg, html=snakemake.output.html)