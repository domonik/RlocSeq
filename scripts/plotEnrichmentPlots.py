import pandas as pd
from DEplots.enrichment import enrichment_plot_from_cp_table
from plotly.subplots import make_subplots
import numpy as np

def main(tsv1, tsv2, tsv3, config, svg, html):
    colorscale = ["purple", "yellow"]
    fig = make_subplots(rows=1, cols=3, shared_yaxes=True, column_titles=["- Puromycin", "+ Puromycin", "intersection"])

    for i, tsv in enumerate((tsv1, tsv2, tsv3)):
        df = pd.read_csv(tsv, sep="\t", )
        df["p.adjust"] = -1 * np.log10(df["p.adjust"])

        f1 = enrichment_plot_from_cp_table(
            df,
            mode="scatter",
            colorscale=colorscale
        )
        if i > 0:
            fig.update_traces(showlegend=False)
        fig.add_traces(f1.data, cols=1+i, rows=1)
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config["EnrichFigure6"]["layout"])
    fig.update_layout(coloraxis=f1.layout.coloraxis)
    fig.update_layout(coloraxis=dict(colorbar=dict(title=dict(text="-log<sub>10</sub>(p.adjust)"))))


    fig.write_image(svg)
    fig.write_html(html)
    p = 0

if __name__ == '__main__':
    tsv1 = snakemake.input.tsv1
    tsv2 = snakemake.input.tsv2
    tsv3 = snakemake.input.tsv3

    main(tsv1, tsv2, tsv3, snakemake.config, snakemake.output.svg, snakemake.output.html)

