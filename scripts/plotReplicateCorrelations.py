import pandas as pd
from DEplots.qc import pheatmap
from scripts.plotly_template import register_template

register_template()

def main(tsv, cfg, layout, outfile):
    df = pd.read_csv(tsv, sep="\t")
    colors = [cfg["plusPColor"], "white", cfg["minusPColor"]]
    fig = pheatmap(df, colorscale=colors, vertical_spacing=0.0, horizontal_spacing=0., row_heights=[0.1, 0.9])
    fig.update_yaxes(showline=False, showgrid=False, zeroline=False, col=2)
    fig.update_yaxes(showline=False, showgrid=False, zeroline=False, row=1, showticklabels=False)
    fig.update_xaxes(showline=False, showgrid=False, zeroline=False, col=2, showticklabels=False)
    fig.update_xaxes(showline=False, showgrid=False, zeroline=False, row=1)

    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(layout)

    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig.update_xaxes(row=1, ticks="", showline=False)
    fig.update_xaxes(row=2, ticks="", showline=False)
    fig.update_yaxes(row=1, ticks="", showline=False)
    fig.update_yaxes(row=2, ticks="", showline=False)
    fig.write_image(outfile)




if __name__ == '__main__':

    if "snakemake" in globals():
        tsv1 = snakemake.input.tsv1
        cfg = snakemake.config
        output = snakemake.output.svg
        layout = snakemake.params['layout']

    else:
        import yaml

        tsv1 = "Pipeline/RUNS/LightRunPlusPuromycin/PipelineData/IntermediateData/CorrData.tsv"


        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)

    main(tsv1, cfg, layout, output)