from helper import compare_results
from scripts.plotly_template import register_template

register_template()

def main(file, params, config, outfile):
    lout = config[params["figureConfig"]]
    fig = compare_results(
        [(params["name"], file)],
        sep="\t",
        condition=params["condition"],
        base=params["baseline"],
        padj_cutoff=params["pAdjCutOff"],
        lfc_cutoff=params["log2FCCutOff"],
        colors=(config["cytosolColor"], "white", config["membraneColor"]),
        y="Name")

    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(lout["layout"])

    #fig.update_layout(config["LayoutDefault"])
    fig.update_yaxes(title=None, showticklabels=False, showline=True, mirror=True, linewidth=1, linecolor="black")
    fig.update_xaxes(title="Fraction of transcriptome [%]", showline=True, mirror=True, linewidth=1,
                     linecolor="black", gridcolor="black", dtick=20)
    fig.update_layout(
        legend=dict(title=None),
        xaxis=dict(range=[-1, 101]),
        yaxis=dict(ticks="")
    )
    fig.write_image(outfile)



if __name__ == '__main__':
    main(
        snakemake.input.tsv,
        snakemake.params,
        snakemake.config,
        snakemake.output.svg
    )