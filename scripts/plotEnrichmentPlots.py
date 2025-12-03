import pandas as pd
from DEplots.enrichment import enrichment_plot_from_cp_table
from plotly.subplots import make_subplots
import numpy as np
from scripts.plotly_template import register_template, PT7

register_template()

def main(tsv1, tsv2, tsv3, config, svg, html, png):
    colorscale = ["purple", "yellow"]
    fig = make_subplots(rows=1, cols=3, shared_yaxes=True, column_titles=["Untreated", "Puromycin", "Intersection"], x_title="Gene ratio")
    fig.update_annotations(font=dict(size=PT7))
    dfs = []
    for i, tsv in enumerate((tsv1, tsv2, tsv3)):
        df = pd.read_csv(tsv, sep="\t", )
        dfs.append(df)
        df["p.adjust"] = -1 * np.log10(df["p.adjust"])

        f1 = enrichment_plot_from_cp_table(
            df,
            mode="scatter",
            colorscale=[config["plusPColor"], config["minusPColor"]]
        )
        if i > 0:
            fig.update_traces(showlegend=False)
        fig.add_traces(f1.data, cols=1+i, rows=1)
    max_ = len(set(dfs[0]["ID"].tolist() + dfs[1]["ID"].tolist() + dfs[2]["ID"].tolist()))
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config["SupplementFullEnrichment"]["layout"])
    fig.update_layout(coloraxis=f1.layout.coloraxis)
    fig.update_yaxes(range=[-0.5, max_ - 0.5])
    fig.update_layout(coloraxis=dict(colorbar=dict(title=dict(text="-log<sub>10</sub>(p.adjust)"))))


    fig.write_image(svg)
    fig.write_image(png, scale=3)
    fig.write_html(html)
    p = 0

if __name__ == '__main__':
    tsv1 = snakemake.input.tsv1
    tsv2 = snakemake.input.tsv2
    tsv3 = snakemake.input.tsv3

    main(tsv1, tsv2, tsv3, snakemake.config, snakemake.output.svg, snakemake.output.html, snakemake.output.png)

