import pandas as pd
from DEplots.enrichment import plot_gsea


def main(df, info, config, output):
    df = pd.read_csv(df, sep="\t")
    info = pd.read_csv(info, sep="\t")
    info.Description = info.Description.str.replace("_", " ").str.replace("and", "&")
    df.Description = df.Description.str.replace("_", " ").str.replace("and", "&")

    fig = plot_gsea(
        df, info, show_zero_lfc=True,
        condition_name="Membrane",
        base_name="Cytosol"
    )
    for d in fig.data[0:2]:
        d.showlegend=False
    fig.add_annotation(
        yref=f"y domain",
        y=0,
        x=1,
        textangle=90,
        xref="x domain",
        xanchor="left",
        yanchor="bottom",
        text=f"p-value:",
        showarrow=False
    )
    for idx, row in info.iterrows():
        pval = row["p.adjust"]
        fig.add_annotation(
            yref=f"y{idx+2} domain",
            y=0.5,
            x=1,
            textangle=90,

            xref="x domain",
            xanchor="left",
            yanchor="middle",
            text=f"{pval:.2e}",
            showarrow=False
        )
    fig.update_shapes(line=dict(width=2))
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config["RBP3GSEA"]["layout"])

    fig.write_image(output.svg)
    fig.write_html(output.html)



if __name__ == '__main__':
    main(
        snakemake.input.df,
        snakemake.input.info,
        snakemake.config,
        snakemake.output
    )