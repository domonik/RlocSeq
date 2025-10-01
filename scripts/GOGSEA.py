from DEplots.enrichment import plot_gsea
import plotly.graph_objects as go
import pandas as pd

def main(
        tsv,
        info_tsv,
        descs,
        layout,
        condition_name: str = None,
        base_name: str = None,

        config:dict = None,
        outfile: str = None,
        **kwargs
):
    df = pd.read_csv(tsv, sep="\t")
    info = pd.read_csv(info_tsv, sep="\t")
    fig: go.Figure = plot_gsea(
        df=df,
        info=info,
        descs=descs,
        show_zero_lfc=True,
        condition_name=condition_name,
        base_name=base_name,
        vertical_spacing=0,


    )
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(layout)
    fig.write_image(outfile)


if __name__ == '__main__':
    main(
        snakemake.input.tsv,
        snakemake.input.info_tsv,
        layout=snakemake.params.layout,
        descs=["ribosome"],
        config=snakemake.config,
        outfile=snakemake.output.svg,
        base_name=snakemake.params.base_name,
        condition_name=snakemake.params.condition_name,
    )