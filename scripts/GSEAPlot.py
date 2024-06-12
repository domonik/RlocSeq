import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
import numpy as np
from plotly.subplots import make_subplots


def plot_gsea(df, info, descs, colors=px.colors.DEFAULT_PLOTLY_COLORS):
    d = descs
    df = df[df["Description"].isin(d)]
    df['Description'] = pd.Categorical(df['Description'], categories=d, ordered=True)
    df = df.sort_values(by=["Description", "x"])
    fig_old = px.line(df, x="x", y="runningScore", color='Description')

    fig = make_subplots(
        rows=len(d) + 1, shared_xaxes=True, vertical_spacing=0,
        row_heights=[0.6] + [0.4 / len(d) for _ in d]

    )
    df[df["position"] == 0]["position"] = np.nan
    for idx, data in enumerate(fig_old.data):
        data.update(marker=dict(color=colors[idx]), line=dict(color=colors[idx]))
        fig.add_trace(data, row=1, col=1)
    for idx, desc in enumerate(d, 1):
        sdf = df[df["Description"] == desc]
        y = sdf["position"]

        fig.add_trace(
            go.Bar(
                x=sdf["x"],
                y=y,
                marker_color=colors[idx-1],
                marker_line=dict(width=1, color=colors[idx-1]),
                name=desc + f' {info[info["Description"] == desc]["p.adjust"].iloc[0]:.2e}',
                hovertext=df["gene"],
            ),
            row=idx + 1, col=1
        )
        fig.update_yaxes(range=[0, 1], row=idx +1, showticklabels=False)
    fig.update_xaxes(title_text="Position in ranked dataset", row=len(d)+1)
    fig.update_yaxes(title_text="Running enrichment score", row=1)

    return fig




if __name__ == '__main__':
    file1 = snakemake.input["df"]
    file2 = snakemake.input["info"]
    df = pd.read_csv(file1, sep="\t")
    descs = sorted(df["Description"].unique())
    info = pd.read_csv(file2, sep="\t")
    fig = plot_gsea(df, info, descs)
    fig.write_html(snakemake.output.html)