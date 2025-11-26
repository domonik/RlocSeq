import pandas as pd
import numpy as np
import plotly.graph_objs as go
import re
from scripts.plotly_template import register_template

register_template()

def extract_terms_before_brackets(text):
    if pd.isna(text):
        return []
    text = text.replace("SUBCELLULAR LOCATION: ", "")
    text = text.split(".")
    ret = [i.split("{")[0].strip() for i in text if i != "" and "Note" not in i]
    print(ret)
    return ret

# Assuming your DataFrame is named 'df' and the column is 'text_column'


def check_membrane_dist(table, config, output, uni_col="Gene Names (ordered locus)"):
    df = pd.read_csv(table, sep="\t")
    x = ("Cytoplasm", "Plasmamembrane", "Cellular thylakoid membrane")

    #odf = odf[["locus_tag_x", "Gene Names", "mRNA", "UTR", "antisense", "non-coding"] + list(x)]

    rdf = df.copy()
    rdf["Transmembrane"] = ~rdf["Transmembrane"].isna()
    rdf = rdf[~rdf["tRNA"]]
    rdf = rdf[rdf["mRNA"]]
    rdf = rdf[["locus_tag", "log2FoldChange_minusPuromycin", "log2FoldChange_plusPuromycin", "Transmembrane", "padj_plusPuromycin", "padj_minusPuromycin"]]
    rdf.to_csv(output.r_is_not_a_valid_programming_language_and_no_one_should_use_this_bullshit, sep="\t", index=False)
    df = df[df["baseMean_minusPuromycin"] > 0]
    df = df[df["mRNA"]]

    df["Membrane enriched"] = df["minusPMembran"]
    df["Cytosol enriched"] = df["minusPCytosol"]
    df["not significant"] = ~(df["minusPCytosol"] | df["minusPMembran"])

    cell_membrane = df[df["Cell membrane"] == True]
    membrane = df[df["Membrane"] == True]
    inner_membrane = df[df["Plasmamembrane"] == True]
    thylakoid = df[df["Cellular thylakoid membrane"] == True]
    cytoplasm = df[df["Cytoplasm"] == True]
    fig = go.Figure()



    all = np.asarray([cytoplasm.shape[0], inner_membrane.shape[0], thylakoid.shape[0]])
    traces = []
    parts = {
        "Cytosol enriched": (config["cytosolColor"], None),
        "not significant": ("rgb(234, 234, 234)", None),
        "Membrane enriched": (config["membraneColor"], None),

    }
    for part, (color, pattern) in parts.items():
        t = np.asarray([cytoplasm[part].sum(), inner_membrane[part].sum(), thylakoid[part].sum()])
        print(t)
        trace = (t / all) * 100
        if np.all(trace == 0):
            continue
        traces.append(trace)
        fig.add_trace(
            go.Bar(
                x=trace,
                y=x,
                name=part,
                orientation="h",
                marker=dict(color=color, pattern_shape=pattern, pattern=dict(fgcolor="grey"), line=dict(color="black", width=1)),
                insidetextanchor="middle",

            )
        )
    fig.update_layout(
        barmode="stack",

    )
    fig.update_layout(
        config["LayoutDefault"]

    )
    fig.update_layout(
        config["F3SubCellLoc"]["layout"]
    )
    fig.update_xaxes(range=[-1, 101])
    fig.show()
    fig.write_image(output.svg, scale=6)



if __name__ == '__main__':
    check_membrane_dist(
        snakemake.input.deseq_res,
        snakemake.config,
        snakemake.output

    )

