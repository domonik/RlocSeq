import pandas as pd
import numpy as np
import plotly.graph_objs as go
import re


def extract_terms_before_brackets(text):
    if pd.isna(text):
        return []
    text = text.replace("SUBCELLULAR LOCATION: ", "")
    text = text.split(".")
    ret = [i.split("{")[0].strip() for i in text if i != "" and "Note" not in i]
    print(ret)
    return ret

# Assuming your DataFrame is named 'df' and the column is 'text_column'


def check_membrane_dist(uniprot_file, deseq_table, config, output, uni_col="Gene Names (ordered locus)"):
    x = ("Cytoplasm", "Membrane", "Cell membrane", "Cell inner membrane", "Cellular thylakoid membrane")

    df = pd.read_csv(uniprot_file, sep="\t")
    subcell = "Subcellular location [CC]"
    df = df.rename({uni_col: "locus_tag"}, axis=1)
    df["locus_tag"] = df["locus_tag"].str.split("; ")
    df = df.explode('locus_tag').reset_index(drop=True)
    print(df)

    df = df[["locus_tag", "Gene Names", subcell, "Transmembrane", "Sequence"]] # dont need intramembrane since there are zero annotated
    df['Cell membrane'] = df[subcell].str.contains("Cell membrane") | df[subcell].str.contains("Cell inner membrane")
    df["Cell inner membrane"] = df[subcell].str.contains("Cell inner membrane")
    df['Cellular thylakoid membrane'] = df[subcell].str.contains("Cellular thylakoid membrane")
    df['Cytoplasm'] = (df[subcell].str.contains("Cytoplasm") == True) & (~(df[subcell].str.contains("Cytoplasmic side") == True))
    df.loc[df[subcell].isna(), "Cytoplasm"] = pd.NA
    df['Membrane'] = (df[subcell].str.contains("Membrane") == True) | df['Cellular thylakoid membrane'] | df['Cell membrane'] | df["Cell inner membrane"]
    df.loc[df[subcell].isna(), "Membrane"] = pd.NA


    deseq_res = pd.read_csv(deseq_table, sep="\t")
    deseq_res = deseq_res.rename({"Unnamed: 0": "locus_tag"}, axis=1)
    deseq_res["merge_tag"] = deseq_res["locus_tag"]
    deseq_res.loc[deseq_res["locus_tag"].str.contains("UTR"), "merge_tag"] = deseq_res.loc[deseq_res["locus_tag"].str.contains("UTR"), "merge_tag"].str[:-5]

    odf = deseq_res.merge(df, left_on="merge_tag", right_on="locus_tag", how="left")
    df = deseq_res.merge(df, left_on="locus_tag", right_on="locus_tag")

    #odf = odf[["locus_tag_x", "Gene Names", "mRNA", "UTR", "antisense", "non-coding"] + list(x)]
    odf = odf.rename({"locus_tag_x": "locus_tag"}, axis=1)
    odf.to_csv(output.tsv, sep="\t", index=False)

    df["-Puromycin Membrane"] = df["minusPMembran"] & ~df["stayedAtMembrane"]
    df["+Puromycin Membrane"] = df["plusPMembran"] & ~df["stayedAtMembrane"]
    df["-Puromycin Cytosol"] = df["minusPCytosol"] & ~df["stayedInCytosol"]
    df["+Puromycin Cytosol"] = df["plusPCytosol"] & ~df["stayedInCytosol"]
    df["nowhere significant"] = ~(df["plusPCytosol"] | df["minusPCytosol"] | df["minusPMembran"] | df["plusPMembran"])

    cell_membrane = df[df["Cell membrane"] == True]
    membrane = df[df["Membrane"] == True]
    inner_membrane = df[df["Cell inner membrane"] == True]
    thylakoid = df[df["Cellular thylakoid membrane"] == True]
    cytoplasm = df[df["Cytoplasm"] == True]
    fig = go.Figure()



    all = np.asarray([cytoplasm.shape[0], membrane.shape[0], cell_membrane.shape[0], inner_membrane.shape[0], thylakoid.shape[0]])
    traces = []
    parts = {
        "stayedInCytosol": ("rgb(0, 38, 149)", None) ,
        "-Puromycin Cytosol": ("rgba(0, 38, 149, 0.8)", "."),
        "+Puromycin Cytosol": ("rgba(0, 38, 149, 0.8)", "x"),
        "nowhere significant": ("rgb(234, 234, 234)", None),
        "+Puromycin Membrane": ("rgba(44, 160, 44, 0.8)", "x"),
        "-Puromycin Membrane": ("rgba(44, 160, 44, 0.8)", "."),
        "stayedAtMembrane": ("rgb(44, 160, 44)", None)

    }
    for part, (color, pattern) in parts.items():
        t = np.asarray([cytoplasm[part].sum(), membrane[part].sum(), cell_membrane[part].sum(), inner_membrane[part].sum(), thylakoid[part].sum()])
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
                marker=dict(color=color, pattern_shape=pattern, pattern=dict(fgcolor="grey"), line=dict(color="black", width=2)),
                text=t,
                textposition="inside",
                insidetextanchor="middle",

            )
        )
    fig.update_layout(
        barmode="stack",

    )
    fig.update_layout(
        **config["CCFigure"]["layout"]
    )
    fig.update_xaxes(range=[-1, 101])

    fig.write_image(output.png, scale=6)



if __name__ == '__main__':
    check_membrane_dist(
        snakemake.input.uniprot,
        snakemake.input.deseq_res,
        snakemake.config,
        snakemake.output

    )

