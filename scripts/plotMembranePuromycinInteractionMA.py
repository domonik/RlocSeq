import numpy as np
import pandas as pd
import plotly.graph_objects as go

from scripts.plotly_template import register_template


register_template()


def main(tsv, gene_names, output, config, params):
    df = pd.read_csv(tsv, sep="\t", index_col=0)
    df.index.name = "locus_tag"
    df = df.reset_index()

    names = pd.read_csv(gene_names)
    df = df.merge(names, how="left", on="locus_tag")
    df["gene_name"] = df["gene_name"].fillna(df["locus_tag"])

    # DESeq2 does not report an interaction estimate for rows without enough
    # information. A logarithmic abundance axis also requires baseMean > 0.
    df = df.loc[
        (df["baseMean"] > 0)
        & np.isfinite(df["baseMean"])
        & np.isfinite(df["log2FoldChange"])
    ].copy()

    fdr_cutoff = config["membrane_puromycin"]["pAdjCutOff"]
    significant = df["padj"].notna() & (df["padj"] <= fdr_cutoff)
    non_significant = df.loc[~significant]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=non_significant["baseMean"],
            y=non_significant["log2FoldChange"],
            mode="markers",
            marker=dict(
                color=config["seahawksGrey"],
                size=3,
                opacity=0.45,
            ),
            customdata=np.stack(
                (
                    non_significant["gene_name"],
                    non_significant["locus_tag"],
                    non_significant["padj"],
                ),
                axis=-1,
            ),
            hovertemplate=(
                "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
                "Mean normalized count: %{x:.2f}<br>"
                "Interaction log2FC: %{y:.2f}<br>"
                "FDR adjusted <i>P</i>: %{customdata[2]:.3g}"
                "<extra></extra>"
            ),
            showlegend=False,
        )
    )

    for is_positive, color, label in (
        (False, config["minusPColor"], "Decreased membrane enrichment"),
        (True, config["plusPColor"], "Increased membrane enrichment"),
    ):
        sdf = df.loc[significant & ((df["log2FoldChange"] > 0) == is_positive)]
        if sdf.empty:
            continue
        fig.add_trace(
            go.Scatter(
                x=sdf["baseMean"],
                y=sdf["log2FoldChange"],
                mode="markers",
                marker=dict(color=color, size=3, opacity=0.9),
                name=label,
                customdata=np.stack(
                    (sdf["gene_name"], sdf["locus_tag"], sdf["padj"]),
                    axis=-1,
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
                    "Mean normalized count: %{x:.2f}<br>"
                    "Interaction log2FC: %{y:.2f}<br>"
                    "FDR adjusted <i>P</i>: %{customdata[2]:.3g}"
                    "<extra></extra>"
                ),
            )
        )

    fig.add_hline(
        y=0,
        line=dict(color="rgb(100, 100, 100)", width=1, dash="dash"),
    )
    fig.add_annotation(
        x=0.98,
        y=0.98,
        xref="paper",
        yref="paper",
        xanchor="right",
        yanchor="top",
        text=f"FDR < {fdr_cutoff:g}: n = {int(significant.sum())}",
        showarrow=False,
    )

    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config[params.figureConfig]["layout"])
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
