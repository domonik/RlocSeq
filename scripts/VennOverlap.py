import pandas as pd
from plotly.subplots import make_subplots
from scripts.plotlyVenn import venn_to_plotly, copy_to_subfigures
import plotly.graph_objs as go


def main(tsv, config, svg):
    df = pd.read_csv(tsv, sep="\t")
    plus_p_m = set(df[df["plusPMembran"]].locus_tag)
    minus_p_m = set(df[df["minusPMembran"]].locus_tag)
    minus_p_c = set(df[df["minusPCytosol"]].locus_tag)
    plus_p_c = set(df[df["plusPCytosol"]].locus_tag)
    multi_fig = make_subplots(rows=1, cols=2, column_titles=["Cytosol", "Membrane"])
    f1 = venn_to_plotly([minus_p_m, plus_p_m, ], L_labels=["Minus Puromycin", "Plus Puromyin"])
    f2 = venn_to_plotly([minus_p_c, plus_p_c,], L_labels=["Minus Puromycin", "Plus Puromyin"])
    figures = [[f2, f1]]
    multi_fig = copy_to_subfigures(figures, cols=2, multi_fig=multi_fig, n_sets=2)

    multi_fig.update_layout(config["LayoutDefault"])
    multi_fig.update_layout(config["Fig5Venn"]["layout"])
    multi_fig.update_xaxes(showgrid=False, zeroline=False)
    multi_fig.update_yaxes(showgrid=False, zeroline=False)
    for annotation in multi_fig.layout.annotations:
        if annotation.text.isnumeric():
            annotation.update(font=dict(size=8))
        else:
            annotation.update(font=dict(size=12))
    multi_fig.write_image(svg)

    p = 0




if __name__ == '__main__':
    file = snakemake.input.tsv


    main(
        file,
        config = snakemake.config,
        svg=snakemake.output.svg
    )