import pandas as pd
from plotly.subplots import make_subplots
from scripts.plotlyVenn import venn_to_plotly
import plotly.graph_objs as go
from scripts.plotly_template import register_template
from scripts.plotFISH import _color_to_calpha

register_template()

def main(file, config, output):
    df = pd.read_csv(file, sep="\t")
    df = df[df["mRNA"]]
    locs = ["All", "Cytoplasm", "Plasmamembrane", "Cellular thylakoid membrane"]
    rows = 2
    cols = 2
    figures = [[] for _ in range(rows)]

    for idx, loc in enumerate(locs):
        if loc == "All":
            sdf = df
        else:
            sdf = df[df[loc] == True]
        rbp3_set = set(sdf[sdf["Rbp3_Target"] == True]["locus_tag"])
        minusp_set = set(sdf[sdf["minusPMembran"] == True]["locus_tag"])
        plusp_set = set(sdf[sdf["plusPMembran"] == True]["locus_tag"])
        subfig = venn_to_plotly(
            (minusp_set, rbp3_set, plusp_set),
            ("Untreated at membrane", "Rbp3 Target", "Puromycin at membrane"),
            L_color=[_color_to_calpha(config["plusPColor"], 0.4), _color_to_calpha(config["minusPColor"], 0.4), _color_to_calpha(config["seahawksGrey"], 0.4)])
        row = idx // cols
        figures[row].append(subfig)
    stitels = ["All", "Cytoplasm", "Plasma membrane", "Thylakoid membrane"]
    multi_fig = make_subplots(rows=rows, cols=cols, horizontal_spacing=0, vertical_spacing=0.05, subplot_titles=stitels)
    multi_fig.update_annotations(font=dict(size=8))
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
            multi_fig.update_xaxes(fig["layout"]["xaxis"], row=row, col=idx)
            multi_fig.update_yaxes(fig["layout"]["yaxis"], row=row, col=idx)
            multi_fig.update_yaxes(scaleanchor="x" if running_idx == 1 else f"x{running_idx}", scaleratio=1, col=idx)
            multi_fig.update_xaxes(scaleanchor="y" if running_idx == 1 else f"x{running_idx}", scaleratio=1, col=idx)
    multi_fig.update_layout(config["LayoutDefault"])
    multi_fig.update_layout(config["RBP3Venn"]["layout"])
    multi_fig.update_xaxes(showgrid=False, zeroline=False, showline=False)
    multi_fig.update_yaxes(showgrid=False, zeroline=False, showline=False)
    for annotation in multi_fig.layout.annotations:
        if annotation.text.isnumeric():
            annotation.update(font=dict(size=8))
        else:
            annotation.update(font=dict(size=10))
    multi_fig.write_html(output.html)
    multi_fig.write_image(output.svg)




if __name__ == '__main__':
    main(
        snakemake.input.tsv,config=snakemake.config, output=snakemake.output
    )
