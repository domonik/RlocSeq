import plotly.graph_objects as go
import yaml
import pandas as pd



def main(tsv, cfg, outfile):
    df = pd.read_csv(tsv, sep="\t")
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=df["Target"],
            y=df["log2FC"],
            error_y=dict(
                type='data',
                visible=True,
                array=df["Std. Dev.1"]
            ),
            marker=dict(line=dict(color="black", width=1), color=cfg["cytosolColor"])
        )
    )
    fig.update_yaxes(title=dict(text="log<sub>2</sub>(fold change)"))
    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(cfg["F4QTPCR"]["layout"])

    fig.write_image(outfile)



if __name__ == '__main__':
    with open("configFile.yml", 'r') as f:
        cfg = yaml.safe_load(f)
    cfg = snakemake.config
    outfile = snakemake.output.svg
    tsv = snakemake.input.tsv
    main(tsv, cfg, outfile)