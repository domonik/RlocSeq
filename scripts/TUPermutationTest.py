import time

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from multiprocessing import Pool
from plotly.subplots import make_subplots
from scripts.plotly_template import register_template

register_template()


def calc_mean_fraction(df, val):
    grouped = df.groupby('shortest_TU')[val].agg(
        fraction=lambda x: x.sum() / len(x),
        count='size'
    ).reset_index()
    grouped = grouped[grouped['count'] > 1]
    original_value = grouped["fraction"].mean()
    return original_value


def mp_wrapper(df):
    perm = np.random.permutation(len(df))
    jdf = df[["minusPMembran", "minusPCytosol"]]

    shuffled = df[['shortest_TU', 'TU_size']].iloc[perm].reset_index(drop=True)
    combined = pd.concat([shuffled.reset_index(drop=True), jdf.reset_index(drop=True)], axis=1)

    value1 = calc_mean_fraction(combined, "minusPMembran")
    value2 = calc_mean_fraction(combined, "minusPCytosol")
    return value1, value2


def main(df, values, cfg,  outfile):
    df = df[df["mRNA"]]
    original_value_mem = calc_mean_fraction(df, 'minusPMembran')
    original_value_cyt = calc_mean_fraction(df, 'minusPCytosol')

    df = df[["minusPMembran", "minusPCytosol", 'shortest_TU', 'TU_size']].copy()

    values = np.load(values)
    s = time.time()
    fig = make_subplots(cols=2, rows=1, subplot_titles=["Membrane", "Cytosol"])
    fig.update_annotations(font=dict(size=8))
    mem_pval = ((values[:, 0] >= original_value_mem).sum() + 1) / (values.shape[0] + 1)
    cyt_pval = ((values[:, 1] >= original_value_cyt).sum() + 1) / (values.shape[0] + 1)
    fig.add_trace(
        go.Histogram(
            x=values[:, 0],
            marker=dict(color=cfg["membraneColor"]),
            histnorm="probability",
            showlegend=False

        ),
        col=1, row=1
    )
    fig.add_trace(
        go.Histogram(
            x=values[:, 1],
            marker=dict(color=cfg["cytosolColor"]),
            histnorm="probability",  # <--- normalize to probability
            showlegend=False

        ),
        col=2, row=1
    )
    fig.add_vline(x=original_value_mem, row=1, col=1, line=dict(color="black", width=1))
    fig.add_vline(x=original_value_cyt, row=1, col=2, line=dict(color="black", width=1))
    fig.add_annotation(
        x=original_value_mem,
        y=1,  # relative y in paper coordinates (1 = top)
        xref="x1",  # first subplot
        yref="paper",
        text=f"{original_value_mem:.3f}",
        showarrow=False,
        yanchor="top",
        xanchor="left"
    )
    fig.add_annotation(
        x=0,
        y=1,  # relative y in paper coordinates (1 = top)
        xref="x domain",  # first subplot
        yref="paper",
        text=f"p-Value: {mem_pval:.1e}",
        showarrow=False,
        yanchor="top",
        xanchor="left",
        bgcolor="white",        # background color of the box
    )
    fig.add_annotation(
        x=0,
        y=1,  # relative y in paper coordinates (1 = top)
        xref="x2 domain",  # first subplot
        yref="paper",
        text=f"p-Value: {cyt_pval:.1e}",
        showarrow=False,
        yanchor="top",
        xanchor="left",
        bgcolor="white",        # background color of the box
    )

    fig.add_annotation(
        x=original_value_cyt,
        y=1,
        xref="x2",  # second subplot
        yref="paper",
        text=f"{original_value_cyt:.3f}",
        showarrow=False,
        yanchor="top",
        xanchor="left",


    )
    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(cfg["SupplementTUAssociation"]["layout"])
    fig.update_layout(xaxis=dict(title="Mean fraction sig. enriched"))
    fig.update_layout(xaxis2=dict(title="Mean fraction sig. enriched"))
    fig.write_image(outfile)
    p = 0



if __name__ == '__main__':
    if "snakemake" in globals():
        tsv = snakemake.input.tsv
        cfg = snakemake.config
        outfile = snakemake.output.svg
        values = snakemake.input.values
    else:
        tsv = "Pipeline/Paper/Supplementary/SummaryTable.tsv"
        import yaml
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        threads = 4
        num_permutations = 100
    df = pd.read_csv(tsv, sep="\t")
    main(df, values=values, outfile=outfile, cfg=cfg)