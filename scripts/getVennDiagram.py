from scripts.plotlyVenn import venn_to_plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import string
import pandas as pd



def plot_venn_diagramm(minus_puro, plus_puro, config, output):

    m_puro = []
    p_puro = []
    for name, file in minus_puro:
        df = pd.read_csv(file, sep="\t")
        df_up = df[(df["log2FoldChange"] >= config["primary"]["log2FCCutOff"]) & (
                df["padj"] <= config["primary"]["pAdjCutOff"])]
        df_down = df[(df["log2FoldChange"] <= -config["primary"]["log2FCCutOff"]) & (
            df["padj"] <= config["primary"]["pAdjCutOff"])]
        m_puro.append((name, df_up, df_down))

    for name, file in plus_puro:
        df = pd.read_csv(file, sep="\t")
        df_up = df[(df["log2FoldChange"] >= config["primary"]["log2FCCutOff"]) & (
                df["padj"] <= config["primary"]["pAdjCutOff"])]
        df_down = df[(df["log2FoldChange"] <= -config["primary"]["log2FCCutOff"]) & (
            df["padj"] <= config["primary"]["pAdjCutOff"])]
        p_puro.append((name, df_up, df_down))
    cols = 2
    rows = 2
    figdown = venn_to_plotly(
        [set(data[2].index) for data in m_puro],
        [data[0] for data in m_puro]
    )
    figup = venn_to_plotly(
        [set(data[1].index) for data in m_puro],
        [data[0] for data in m_puro]

    )
    figupplus = venn_to_plotly(
        [set(data[1].index) for data in p_puro],
        [data[0] for data in p_puro]

    )
    figdownplus = venn_to_plotly(
        [set(data[2].index) for data in p_puro],
        [data[0] for data in p_puro]

    )
    multi_fig = make_subplots(rows=rows, cols=cols, column_titles=("Membrane", "Cytosol"), row_titles=("minus Puromycin", "plus Puromycin"), horizontal_spacing=0, vertical_spacing= 0)
    figures = [(figup, figdown), (figupplus, figdownplus)]
    annos = string.ascii_uppercase

    for row, subfigs in enumerate(figures, 1):
        for idx, fig in enumerate(subfigs, 1):
            running_idx = (row -1) * cols + idx
            for iidx, shape in enumerate(fig['layout']['shapes']):
                shape["xref"] = f"x{running_idx}" if running_idx > 1 else "x"
                shape["yref"] = f"y{running_idx}" if running_idx > 1 else "y"
                multi_fig.add_shape(shape)
                color = shape["fillcolor"]
                if running_idx == 1 and color != "rgba(255,0,0,0)":
                    multi_fig.add_trace(
                        go.Scatter(
                            x=[None],
                            y=[None],
                            mode="markers",
                            name=fig['layout']['annotations'][iidx]["text"],
                            marker=dict(
                                color=color,
                                size=900,
                                # line=dict(
                                #     color="black",
                                #     width=3
                                # ),
                            ),
                            showlegend=True

                        )
                    )
            for annotation in fig['layout']['annotations'][3:]:
                if annotation["xref"] == "x":
                    annotation["xref"] = f"x{running_idx}" if running_idx > 1 else "x"
                    annotation["yref"] = f"y{running_idx}" if running_idx > 1 else "y"
                else:
                    annotation["xref"] = f"x{running_idx} domain" if running_idx > 1 else "x domain"
                    annotation["yref"] = f"y{running_idx} domain" if running_idx > 1 else "y domain"
                multi_fig.add_annotation(annotation)
            multi_fig.add_annotation(
                text=f"<b>{annos[running_idx - 1]}</b>",
                xref=f"x{running_idx} domain" if running_idx > 1 else "x domain",
                xanchor="left",
                x=0,
                yref=f"y{running_idx} domain" if running_idx > 1 else "y domain",
                yanchor="top",
                y=1,
                font=dict(size=config["multipanel_font_size"]),
                showarrow=False

            )
            multi_fig.update_xaxes(fig["layout"]["xaxis"], row=row, col=idx)
            multi_fig.update_yaxes(fig["layout"]["yaxis"], row=row, col=idx)
            multi_fig.update_yaxes(scaleanchor="x" if running_idx == 1 else f"x{running_idx}", scaleratio=1, col=idx)
            multi_fig.update_xaxes(scaleanchor="y" if running_idx == 1 else f"x{running_idx}", scaleratio=1, col=idx)
    fig.write_image(output.svg)
    fig.write_html(output.html)


if __name__ == '__main__':
    plot_venn_diagramm(
        [
            ("Light", snakemake.input.deseq1),
            ("Light-Dark, Light", snakemake.input.deseq2),
            ("Light-Dark, Dark", snakemake.input.deseq3),
        ],
        [
            ("Light", snakemake.input.lp),
            ("Light-Dark, Light", snakemake.input.ldlp),
            ("Light-Dark, Dark", snakemake.input.lddp),
        ],
        snakemake.config,
        snakemake.output
    )
