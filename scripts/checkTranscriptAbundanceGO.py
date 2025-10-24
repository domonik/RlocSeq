import pandas as pd
import plotly.graph_objects as go



def main(tsv, godb):
    df = pd.read_csv(tsv, sep="\t")
    godb = pd.read_csv(godb, sep="\t")
    gids = godb[godb.Description == "plasma membrane"].iloc[0].universeGeneID.split("/")
    fig = go.Figure()
    lines = df[df.index.isin(gids)]
    fig.add_trace(
        go.Violin(
            x=["bla" for _ in range(len(lines))],
            y=lines["log2FoldChange"],
        )
    )
    breakpoint()



if __name__ == '__main__':
    main(snakemake.input.tsv, snakemake.input.go)