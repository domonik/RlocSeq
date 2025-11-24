import plotly.graph_objects as go
import yaml
import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
from scripts.plotly_template import register_template

register_template()

def significance_label(p):
    if pd.isna(p): return ""
    elif p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return ""

def main(tsv, cfg, outfile):
    df = pd.read_csv(tsv, sep="\t")
    reference = "16S rRNA"
    ref_values = df.loc[df["Probe"] == reference, "log2FoldChange"].dropna()
    fig = go.Figure()
    order = ["16S rRNA", "rnpB", "psaC", "rimO", "psbA2", "rps21", "rpl10", "rpl17", "psaA", "pilA9-11", "fusB", "cupA"]
    df = df[df["Probe"].isin(order)]
    df["Probe"] = pd.Categorical(df["Probe"], categories=order, ordered=True)
    df = df.sort_values("Probe")

    results = []
    for probe in df["Probe"].unique():
        if probe == reference:
            continue
        probe_values = df.loc[df["Probe"] == probe, "log2FoldChange"].dropna()
        if len(probe_values) == 0:
            continue

        t_stat, p_val = ttest_ind(probe_values, ref_values, equal_var=False)
        results.append({"Probe": probe, "t_stat": t_stat, "p_value": p_val})
    ttest_df = pd.DataFrame(results)
    agg = (
        df.groupby("Probe", as_index=False)
        .agg(mean_log2FC=("log2FoldChange", "mean"),
             std_log2FC=("log2FoldChange", "std"))
    )
    agg = agg.merge(ttest_df, on="Probe", how="left")
    agg["significance"] = agg["p_value"].apply(significance_label)
    agg["Probe"] = "<i>" + agg["Probe"] + "</i>"
    fig.add_trace(
        go.Bar(
            x=agg["Probe"],
            y=agg["mean_log2FC"],
            error_y=dict(
                type="data",
                array=agg["std_log2FC"],
                visible=True,
                thickness=1,  # line width of the error bars (in px)
            ),
            marker=dict(line=dict(color="black", width=1),
                        color="white"),
            showlegend=False# show significance above bars
        )
    )
    y_max = (agg["mean_log2FC"] + agg["std_log2FC"]).max()
    offset = y_max * 0.1
    fig.add_trace(
        go.Scatter(
            x=agg["Probe"],
            y=np.full(len(agg), y_max + offset),  # same height for all
            text=agg["significance"],
            mode="text",
            textfont=dict(size=9.33),
            showlegend=False
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