import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
import numpy as np
from plotly.subplots import make_subplots



def plot_gsea(dfs, infos, colors=px.colors.qualitative.Light24):
    d = ["plasma membrane"]
    fig = make_subplots(
        rows=len(d) + 1, shared_xaxes=True, vertical_spacing=0,
        row_heights=[0.6] + [0.4 / len(d) for _ in d]

    )
    for dfidx, df in enumerate(dfs):
        info = infos[dfidx]
        df = df[df["Description"].isin(d)]
        df['Description'] = pd.Categorical(df['Description'], categories=d, ordered=True)
        df = df.sort_values(by=["Description", "x"])
        fig_old = px.line(df, x="x", y="runningScore", color='Description')


        df[df["position"] == 0]["position"] = np.nan
        for idx, data in enumerate(fig_old.data):
            data.update(marker=dict(color=colors[dfidx]), line=dict(color=colors[dfidx]))
            fig.add_trace(data, row=1, col=1)
        for idx, desc in enumerate(d, 1):
            sdf = df[df["Description"] == desc]
            y = sdf["position"]

            fig.add_trace(
                go.Bar(
                    x=sdf["x"],
                    y=y,
                    marker_color=colors[dfidx],
                    marker_line=dict(width=1, color=colors[dfidx]),
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
    enriched1 = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/LightRunPlusPuromycin/PipelineData/Enrichment/GSEAGO_cM_vs_bC.tsv"
    enriched2 = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/LightRunMinusPuromycin/PipelineData/Enrichment/GSEAGO_cM_vs_bC.tsv"
    gsdata1 = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/LightRunPlusPuromycin/PipelineData/Enrichment/GSEAGO_plot_data_cM_vs_bC.tsv"
    gsdata2 = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/LightRunMinusPuromycin/PipelineData/Enrichment/GSEAGO_plot_data_cM_vs_bC.tsv"
    infos = [pd.read_csv(file1, sep="\t") for file1 in [enriched1, enriched2]]
    dfs = [pd.read_csv(file1, sep="\t") for file1 in [gsdata1, gsdata2]]
    fig = plot_gsea(dfs, infos)
    fig.write_html("foo.html")