import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import re

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
    print(text)
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




def main(motifs_file, anno_file, fastas, outfile):
    motifs = pd.read_csv(motifs_file, sep="\t", comment="#")
    anno = pd.read_csv(anno_file, sep="\t", comment="#")
    motifs = motifs.sort_values(by="motif_ALT_ID")

    seq = {
        "locus_tag": [],
        "Nt_sequence": []
    }
    for fasta in fastas:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seq["locus_tag"].append(seq_record.id)
            seq["Nt_sequence"].append(str(seq_record.seq))

    seq = pd.DataFrame(seq)
    df = motifs.merge(anno, left_on="seq_ID", right_on="locus_tag")

    # df = df.merge(seq, left_on="locus_tag", right_on="locus_tag")
    # df["ProteinSeqFromRNA"] = df["Nt_sequence"].apply(translate)
    #
    # test = df["ProteinSeqFromRNA"] == df["Sequence"]
    # print(df[~test][["ProteinSeqFromRNA", "Sequence"]])
    # df[~test].loc["ProteinSeqFromRNA"] = pd.NA
    # df = df.apply(extract_protein_site, axis=1)
    # df["Transmembrane_intervals"] = df["Transmembrane"].apply(extract_transmem_ranges)
    # df["overlap"] = df.apply(check_if_hit_inside_membrane, axis=1)
    #
    # df["UFreq"] = df["Nt_sequence"].apply(lambda x: x.count("U") / len(x))
    df.to_csv(outfile, sep="\t", index=False)





if __name__ == '__main__':
    import os
    file = os.path.join(snakemake.input.dir, "sequences.tsv")
    main(
        motifs_file=file,
        anno_file=snakemake.input.tsv[0],
        fastas=[snakemake.input.fasta, snakemake.input.bg_fasta],
        outfile=snakemake.output.tsv
    )