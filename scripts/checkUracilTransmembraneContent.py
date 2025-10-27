import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import re
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import plotly.io as pio


pio.renderers.default = "browser"

def translate(seq):
    rna = Seq(seq.replace("T", "U"))
    return str(rna.translate())[:-1]


def extract_protein_site(row):
    start = row["start"] - 1
    seq = row["Nt_sequence"]
    assert seq[start] == row["matched_sequence"][0], f'{seq[start:start+5]}, {row["matched_sequence"][0:5]} unequal'
    protein_start = start / 3
    protein_end = row["stop"] / 3
    row["ProteinStart"] = protein_start

    row["ProteinEnd"] = protein_end
    seq = row["ProteinSeqFromRNA"][int(np.floor(protein_start)):int(np.ceil(protein_end))]
    print(seq)
    row["ProteinSeq"] = seq

    return row

def extract_transmem_ranges(text):
    # Find all occurrences of the pattern
    if pd.isna(text):
        return pd.NA
    matches = re.findall(r'TRANSMEM (\d+)\.\.(\d+)', text)
    # Convert them into a list of tuples with integers
    return [(int(start), int(end)) for start, end in matches]

def intervals_overlap(interval1, interval2):
    (a, b) = interval1
    (c, d) = interval2
    return a < d and c < b

def check_if_hit_inside_membrane(row):
    intervals = row["Transmembrane_intervals"]
    if not isinstance(intervals, list):
        return pd.NA
    start = row["ProteinStart"]
    end = row["ProteinEnd"]
    overlap = False
    for interval in intervals:
        if intervals_overlap(interval, (start, end)):
            overlap = True
    return overlap

def get_rna_indices(intervals):
    if not isinstance(intervals, list):
        return pd.NA
    rna_intervals = []
    for interval in intervals:
        start, end = interval
        rna_start = (start - 1) * 3
        rna_end = (end - 1) * 3
        rna_intervals.append([rna_start, rna_end])
    return rna_intervals

def uracil_content(row):
    seq = row["RNASequence"]
    return seq.count("U") / len(seq)

def test_with_offset(row):
    return row[""]

def main(tsv, outputs, config):
    df = pd.read_csv(tsv, sep="\t", comment="#")
    df = df[~pd.isna(df["Transmembrane"])]
    df = df[df["mRNA"]]


    df["ProteinSeqFromRNA"] = df["RNASequence"].apply(translate)

    test = df["ProteinSeqFromRNA"].str[1:] == df["Sequence"].str[1:]
    df = df[test]
    df["Transmembrane_intervals"] = df["Transmembrane"].apply(extract_transmem_ranges)
    df["RNA_intervals"] = df["Transmembrane_intervals"].apply(get_rna_indices)
    df["Uracil_content"] = df.apply(uracil_content, axis=1)

    #df = df[~pd.isna(df["Uracil_content"])]
    treats = ["minusPuromycin", "plusPuromycin"]
    fig = make_subplots(rows=1, cols=2, column_titles=["- Puromycin", "+ Puromycin"], shared_yaxes=True, x_title="log2FoldChange")
    for idx, puro in enumerate(treats):
        filter = "minusP" if puro == "minusPuromycin" else "plusP"

        for location in ["Cytosol", "Membran"]:
            sdf = df[df[f"{filter}{location}"]]
            color = config["membraneColor"] if location == "Membran" else config["cytosolColor"]
            fig.add_trace(
                go.Scatter(
                    x=sdf[f"log2FoldChange_{puro}"],
                    y=sdf["Uracil_content"],
                    mode="markers",
                    marker=dict(color=color),
                    name=location,
                    customdata=sdf[["gene_name", f"baseMean_{puro}"]],
                    hovertemplate=(
                        "Gene: %{customdata[0]}<br>"
                        "log2FoldChange: %{x:.2f}<br>"
                        "Uracil_content: %{y:.2f}<br>"
                        f"baseMean_{puro}:" + "%{customdata[0]:.2f}<extra>{location}</extra>"
                    ),
                    showlegend=bool(idx)
                ), col=idx+1, row=1

            )
        sdf = df[~(df[f"{filter}Membran"] | df[f"{filter}Cytosol"])]
        fig.add_trace(
            go.Scatter(
                x=sdf[f"log2FoldChange_{puro}"],
                y=sdf["Uracil_content"],
                mode="markers",
                marker=dict(color="#8A8D8F"),
                name="not sig.",
                customdata=sdf[["gene_name", f"baseMean_{puro}"]],
                showlegend=bool(idx),
                hovertemplate=(
                        "Gene: %{customdata[0]}<br>"
                        "log2FoldChange: %{x:.2f}<br>"
                        "Uracil_content: %{y:.2f}<br>"
                        f"baseMean_{puro}:" + "%{customdata[0]:.2f}<extra>{location}</extra>"
                )
            ), col=idx + 1, row=1

        )
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(config["SupplementUracilTransmembrane"]["layout"])
    fig.update_layout(

        yaxis=dict(title="Uracil content [%]")
    )
    fig.write_image(outputs.svg)
    fig.write_html(outputs.html)






if __name__ == '__main__':

    main(
        tsv=snakemake.input.tsv,
        outputs=snakemake.output,
        config=snakemake.config
    )