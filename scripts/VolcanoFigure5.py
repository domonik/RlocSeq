import pandas as pd
from DEplots.volcano import volcano_from_deseq_result
import plotly.graph_objs as go
import numpy as np

def main(tsv, add_data, go_terms, output, config, params):
    df = pd.read_csv(tsv, sep="\t")
    add_data = pd.read_csv(add_data, sep=",")
    df = df.merge(add_data, left_index=True, right_on="locus_tag")
    df = df.set_index("locus_tag")
    go_df = pd.read_csv(go_terms, sep="\t")
    highlight = {}
    vc = params["figureConfig"]
    for term, color in config[vc]["highlight_go_term"].items():
        up_ids = go_df[go_df["Description"] == term]["geneID"].str.split("/").iloc[0]
        universe_ids = go_df[go_df["Description"] == term]["universeGeneID"].str.split("/").iloc[0]
        ids = universe_ids + up_ids
        highlight[term] = (color, ids)
    lout = config[vc]
    fig = volcano_from_deseq_result(
        df,
        name_col="gene_name",
        highlight=highlight,
        lfc_cutoff=config["RUNS"]["secondary"]["log2FCCutOff"],
        padj_cutoff=config["RUNS"]["secondary"]["pAdjCutOff"],
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
        snakemake.input.add_data,
        snakemake.input.go,
        snakemake.output,
        snakemake.config,
        snakemake.params
    )