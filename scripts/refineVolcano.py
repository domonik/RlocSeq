import pandas as pd
from DEplots.volcano import volcano_from_deseq_result
import plotly.graph_objs as go
import numpy as np


def main(tsv, gsea, output, params, config):
    lout = config[params["figureConfig"]]
    df = pd.read_csv(tsv, sep="\t")
    gsea = pd.read_csv(gsea, sep="\t")
    highlight = {key: (color, df[df[key] == True].index) for key, color in config["RNA_TYPE_COLORS"].items() if key in gsea.ID.tolist()}
    df = df.rename({
        "log2FoldChange_minusPuromycin": "log2FoldChange",
        "padj_minusPuromycin": "padj"
    },
    axis=1)
    fig = volcano_from_deseq_result(
        df,
        name_col="gene_name",
        highlight=highlight,
        lfc_cutoff=config["primary"]["log2FCCutOff"],
        padj_cutoff=config["primary"]["pAdjCutOff"],
        condition_name="Membrane",
        base_name="Cytosol",
        highlight_up_color=config["membraneColor"],
        highlight_down_color=config["cytosolColor"],
        opacity=0.3
    )
    fig: go.Figure
    if "annotations" in lout:
        for anno, data in lout["annotations"].items():
            sdf = df[df["gene_name"] == anno]
            print(data)
            for _, n in sdf.iterrows():
                fig.add_annotation(
                    x=n["log2FoldChange"],
                    y=-1 * np.log10(n["padj"]),
                    ax=data["x"],
                    ay=data["y"],
                    text=anno,
                    showarrow=True,
                    axref="x",
                    ayref="y",
                    arrowwidth=1,
                    arrowsize=0.5,
                    arrowcolor="black",
                    font=dict(size=lout["annotationsize"])

                )
    fig.update_layout(config["LayoutDefault"])
    fig.update_xaxes(dtick=2)

    fig.update_layout(lout["layout"])
    fig.update_traces(
        marker=dict(size=3),
    )
    fig.data[0].update(name="not sig.", showlegend=True)
    fig.write_image(output.svg)
    fig.write_html(output.html)

if __name__ == '__main__':
    main(
        snakemake.input.tsv,
        snakemake.input.gsea,
        snakemake.output,
        snakemake.params,
        snakemake.config
    )