import pandas as pd
import numpy as np
import plotly.graph_objs as go


def check_membrane_dist(uniprot_file, deseq_table, config, x=("down", "Not significant", "up"), uni_col="Gene Names (ordered locus)"):
    df = pd.read_csv(uniprot_file, sep="\t")
    deseq_res = pd.read_csv(deseq_table, sep="\t")
    df = deseq_res.merge(df, left_index=True, right_on=uni_col)
    for e in ("Transmembrane", "Intramembrane"):
        df[e] = ~df[e].isna()
    sig = df[df["padj"] <= config["padj_cutoff"]]
    rest = df[df["padj"] > config["padj_cutoff"]]
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]
    fig = go.Figure()
    print(df)
    all = np.asarray([down.shape[0], rest.shape[0], up.shape[0]])
    traces = []
    for part in ("Transmembrane", "Intramembrane"):
        trace = (np.asarray([down[part].sum(), rest[part].sum(), up[part].sum()]) / all) * 100
        traces.append(trace)
        fig.add_trace(
            go.Bar(
                x=x,
                y=trace,
                name=part
            )
        )
    all = np.ones(3) * 100 - np.asarray(traces).sum(axis=0)
    fig.add_trace(
        go.Bar(
            x=x,
            y=all,
            name="Remaining"
        )
    )
    fig.update_layout(
        barmode="stack"
    )
    fig.show()



if __name__ == '__main__':
    check_membrane_dist(
        snakemake.input.uniprot,
        snakemake.input.deseq_res,
        snakemake.config

    )

