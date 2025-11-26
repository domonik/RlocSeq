import pandas as pd
from plotly.subplots import make_subplots
from scripts.plotlyVenn import venn_to_plotly, copy_to_subfigures
import plotly.graph_objs as go
from scipy.stats import pearsonr
from scripts.plotly_template import register_template
from scripts.plotFISH import _color_to_calpha

register_template()

def main(tsv, config, svg):
    df = pd.read_csv(tsv, sep="\t")
    plus_p_m = set(df[df["plusPMembran"]].locus_tag)
    minus_p_m = set(df[df["minusPMembran"]].locus_tag)
    minus_p_c = set(df[df["minusPCytosol"]].locus_tag)
    plus_p_c = set(df[df["plusPCytosol"]].locus_tag)
    multi_fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=["Cytosol", None, "Membrane", None],
        column_widths=[0.3, 0.7],
        vertical_spacing=0.,
        horizontal_spacing=0.2,
    specs=[
        [{}, {"rowspan": 2}],           # First row: normal 2 columns
        [{}, None]  # Second row: one subplot spanning both columns
    ]
    )
    for annotation in multi_fig.layout.annotations:
        annotation.update(font=dict(size=8, ), yanchor="middle")

    colors = [_color_to_calpha(config["minusPColor"], 0.6), _color_to_calpha(config["plusPColor"], 0.6)]
    f1 = venn_to_plotly([minus_p_m, plus_p_m, ], L_labels=["untreated", "Puromycin"], L_color=colors)
    f2 = venn_to_plotly([minus_p_c, plus_p_c,], L_labels=["untreated", "Puromycin"], L_color=colors)
    figures = [[f2,], [f1]]
    multi_fig = copy_to_subfigures(figures, cols=2, multi_fig=multi_fig, n_sets=2)

    multi_fig.update_layout(config["LayoutDefault"])
    multi_fig.update_layout(config["Fig5Venn"]["layout"])
    multi_fig.update_xaxes(showgrid=False, zeroline=False)
    multi_fig.update_yaxes(showgrid=False, zeroline=False)
    white_annos = ["80", "512"]
    for annotation in multi_fig.layout.annotations:
        if annotation.text.isnumeric():
            annotation.update(font=dict(size=8))
            if annotation.text in white_annos:
                annotation.update(font=dict(color="white"))
            elif annotation.text == "135":
                annotation.update(showarrow=True, ax=annotation.x+7, ay=annotation.y-25, arrowcolor="black", arrowwidth=.5, arrowhead=0)
    multi_fig.add_trace(
        go.Scatter(
            x=df["log2FoldChange_minusPuromycin"],
            y=df["log2FoldChange_plusPuromycin"],
            mode="markers",
            showlegend=False,
            marker=dict(size=3, color=cfg["seahawksGrey"])

        ),
        row=1, col=2
    )
    for shape in multi_fig.layout.shapes:
        if shape.xref == "x":
            shape.update(xref="x1", yref="y1")
    multi_fig.update_xaxes(
        title=dict(text="untreated<br>(log<sub>2</sub>FoldChange)"),
        scaleanchor=None,
        range=[-5, 3],
        dtick=1,
        col=2,
        showline=True,

    )
    multi_fig.update_xaxes(
        showline=False,
        col=1
    )
    multi_fig.update_yaxes(
        showline=False,
        col=1
    )
    multi_fig.update_yaxes(
        title=dict(text="+ Puromycin<br>(log<sub>2</sub>FoldChange)"),
        scaleanchor=None,
        col=2,
        dtick=1,
        showline=True,


    )
    mask = ~(df["log2FoldChange_plusPuromycin"].isna() | df["log2FoldChange_minusPuromycin"].isna())
    stat, p = pearsonr(df.loc[mask, "log2FoldChange_minusPuromycin"], df.loc[mask, "log2FoldChange_plusPuromycin"])
    multi_fig.add_annotation(
        text=f"Pearson R: {stat:.2f}",
        xref="x2 domain",
        yref="y2 domain",
        xanchor="right",
        yanchor="bottom",
        x=1,
        y=0,
        showarrow=False,
        bgcolor="white"

    )
    #multi_fig.show()
    multi_fig.write_image(svg)

    p = 0




if __name__ == '__main__':

    if "snakemake" in globals():
        file = snakemake.input.tsv
        cfg = snakemake.config
        outfile = snakemake.output.svg

    else:
        import yaml
        import plotly.io as pio
        pio.renderers.default = "browser"
        file = "Pipeline/Summary/joinedDESeq_cM_vs_bC.tsv"
        outfile = "foo.svg"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)


    main(
        file,
        config = cfg,
        svg=outfile
    )