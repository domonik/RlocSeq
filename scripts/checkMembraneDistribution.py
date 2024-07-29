import pandas as pd
import numpy as np
import plotly.graph_objs as go


def check_membrane_dist(uniprot_file, deseq_table, config, outtable, x=("stayedInCytosol", "Rest", "stayedAtMembrane"), uni_col="Gene Names (ordered locus)"):
    df = pd.read_csv(uniprot_file, sep="\t")
    deseq_res = pd.read_csv(deseq_table, sep="\t")
    deseq_res = deseq_res.rename({"Unnamed: 0": "locus_tag"}, axis=1)

    deseq_res = deseq_res[~deseq_res.locus_tag.str.contains("UTR")]
    deseq_res = deseq_res[~deseq_res.locus_tag.str.contains("-as")]
    deseq_res = deseq_res[~deseq_res.locus_tag.str.contains("ncr|ncl")]
    print(deseq_res)
    df = deseq_res.merge(df, left_on="locus_tag", right_on=uni_col)
    for e in ("Transmembrane", "Intramembrane"):
        df[e] = ~df[e].isna()
    odf = df[["locus_tag", "Transmembrane", "Intramembrane"]]
    odf.to_csv(outtable, sep="\t", index=False)
    up = df[df["stayedAtMembrane"]]
    down = df[df["stayedInCytosol"]]
    rest = df[(~df["stayedAtMembrane"]) & (~df["stayedInCytosol"])]
    fig = go.Figure()
    print(up)
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
        snakemake.config,
        snakemake.output.membrane_table

    )

