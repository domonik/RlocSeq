import pandas as pd
from DEplots.enrichment import plot_gsea
import plotly.express as px
from scripts.plotly_template import register_template

register_template()


def main(df, info, config, output, params):
    lout = config[params["figureConfig"]]
    df = pd.read_csv(df, sep="\t")
    info = pd.read_csv(info, sep="\t")
    info.Description = info.Description.str.replace("_", " ").str.replace("and", "&")
    df.Description = df.Description.str.replace("_", " ").str.replace("and", "&")
    descs = df["Description"].unique()
    to_display = len(descs)

    rh = [0.7] + [0.3 / to_display for _ in range(to_display)]
    colors = [config["cytosolColor"], config["membraneColor"],  "#002244", "black", "red"] if "tRNA" not in descs else [config["RNA_TYPE_COLORS"][key] for key in descs]
    fig = plot_gsea(
        df, info, show_zero_lfc=params["show_zero_lfc"],
        condition_name=params["condition_name"],
        base_name=params["base_name"],
        row_heights=rh,
        vertical_spacing=0,
        colors=colors
    )
    fig.add_annotation(
        yref=f"y domain",
        y=0,
        x=1,
        textangle=90,
        xref="x domain",
        xanchor="left",
        yanchor="bottom",
        text=f"FDF adjusted <i>P</i>:",
        showarrow=False
    )
    for idx, row in info.iterrows():
        pval = row["p.adjust"]
        fig.add_annotation(
            yref=f"y{idx+2} domain",
            y=0.5,
            x=1,
            textangle=0,

            xref="x domain",
            xanchor="left",
            yanchor="middle",
            text=f"{pval:.2e}",
            showarrow=False
        )
    fig.update_shapes(line=dict(width=2))
    fig.update_traces(marker_line=dict(width=0), width=1, selector=dict(type='bar'))
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(lout["layout"])
    fig.layout.annotations = [annotation for annotation in fig.layout.annotations if annotation.text != "0 L2FC"]
    fig.update_traces(
        showlegend=False,  # hide legend for these traces
    )
    fig.update_traces(
        selector=dict(type="bar"),  # only line traces
        showlegend=True,  # hide legend for these traces
        width=2
    )
    fig.update_layout(bargap=0., barmode='overlay')
    for i in range(1, len(rh)):
        fig.update_xaxes(row = i, ticks="")
        fig.update_yaxes(row = i+1, ticks="")
    fig.write_image(output.svg)
    fig.write_html(output.html)



if __name__ == '__main__':
    main(
        snakemake.input.df,
        snakemake.input.info,
        snakemake.config,
        snakemake.output,
        snakemake.params
    )