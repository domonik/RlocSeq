import pandas as pd
import plotly.graph_objs as go

from DEplots.volcano import volcano_from_deseq_result
from scripts.plotly_template import register_template


register_template()


def main(tsv, gene_names, output, config, params):
    df = pd.read_csv(tsv, sep="\t", index_col=0)
    df.index.name = "locus_tag"
    df = df.reset_index()

    names = pd.read_csv(gene_names)
    df = df.merge(names, how="left", on="locus_tag")
    df["gene_name"] = df["gene_name"].fillna(df["locus_tag"])

    fig = volcano_from_deseq_result(
        df,
        name_col="gene_name",
        lfc_cutoff=config["membrane_puromycin"]["log2FCCutOff"],
        padj_cutoff=config["membrane_puromycin"]["pAdjCutOff"],
        condition_name=params.condition_name,
        base_name=params.baseline_name,
        highlight_up_color=config["plusPColor"],
        highlight_down_color=config["minusPColor"],
        opacity=0.05,
    )

    fig: go.Figure
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config[params.figureConfig]["layout"])
    fig.update_annotations(font=dict(color="black"))
    fig.update_xaxes(dtick=2)
    fig.update_traces(marker=dict(size=3))

    for shape in fig.layout.shapes:
        if shape.type == "line":
            shape.update(line=dict(width=1))
        else:
            shape.update(line=dict(width=0))

    fig.data[0].update(name="not sig.", showlegend=False)
    fig.write_image(output.svg)
    fig.write_html(output.html)


if __name__ == "__main__":
    main(
        snakemake.input.tsv,
        snakemake.input.gene_names,
        snakemake.output,
        snakemake.config,
        snakemake.params,
    )
