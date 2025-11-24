import pandas as pd
from DEplots.volcano import volcano_from_deseq_result
import plotly.graph_objs as go
import numpy as np

from scripts.plotly_template import register_template

register_template()

def main(tsv, go_terms, output, config, params):
    suffix = params.suffix
    df = pd.read_csv(tsv, sep="\t")
    df = df[["locus_tag", "gene_name", f"log2FoldChange_{suffix}", f"padj_{suffix}"]]
    df = df.rename({
        f"log2FoldChange_{suffix}": "log2FoldChange",
        f"padj_{suffix}": "padj",
    }, axis=1)
    go_df = pd.read_csv(go_terms, sep="\t")
    highlight = {}
    vc = params["figureConfig"]
    for term, color in config[vc]["highlight_go_term"].items():
        up_ids = go_df[go_df["Description"] == term]["geneID"].str.split("/").iloc[0]
        up_ids = df[df["locus_tag"].isin(up_ids)].index
        universe_ids = go_df[go_df["Description"] == term]["universeGeneID"].str.split("/").iloc[0]
        universe_ids = df[df["locus_tag"].isin(universe_ids)].index

        ids = universe_ids.tolist() + up_ids.tolist()
        highlight[term] = (color, ids)
    lout = config[vc]
    fig = volcano_from_deseq_result(
        df,
        name_col="gene_name",
        highlight=highlight,
        lfc_cutoff=config["RUNS"]["secondary"]["log2FCCutOff"],
        padj_cutoff=config["RUNS"]["secondary"]["pAdjCutOff"],
        condition_name=params.condition_name,
        base_name=params.basename,
        highlight_up_color=config["membraneColor"],
        highlight_down_color=config["cytosolColor"],
        opacity=0.15
    )
    fig.update_annotations(font=dict(color="black"))
    for shape in fig.layout.shapes:
        if shape.type == "line":
            shape.update(line=dict(width=1))
        else:
            shape.update(line=dict(width=0))
    fig: go.Figure
    if "annotations" in lout:
        for anno, data in lout["annotations"].items():
            sdf = df[df["gene_name"] == anno]
            print(data)
            for _, n in sdf.iterrows():
                if "UTR" not in n.locus_tag:
                    fig.add_annotation(
                        x=n["log2FoldChange"],
                        y=-1 * np.log10(n["padj"]),
                        ax=data["x"],
                        ay=data["y"],
                        text=f"<i>{anno}</i>",
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
    fig.data[0].update(name="not sig.", showlegend=False)
    fig.write_image(output.svg)
    fig.write_html(output.html)

if __name__ == '__main__':
    main(
        snakemake.input.tsv,
        snakemake.input.go,
        snakemake.output,
        snakemake.config,
        snakemake.params
    )