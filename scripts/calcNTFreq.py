import pandas as pd
from Bio import SeqIO
import plotly.express as px
from collections import defaultdict
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

NTS = ["A", "U", "C", "G"]

def calc_freqs(seq):
    p = []
    for nt in NTS:
        p.append(seq.count(nt) / len(seq))
    return p


from collections import defaultdict


def trinucleotide_frequency(sequence):
    trinuc_count = defaultdict(int)
    total_trinucs = 0

    # Count trinucleotides
    for i in range(0, len(sequence)):  # Subtract 2 to prevent out-of-range
        trinuc = sequence[i]
        trinuc_count[trinuc] += 1
        total_trinucs += 1

    # Calculate frequencies
    trinuc_freq = {trinuc: count / total_trinucs for trinuc, count in trinuc_count.items()}


    return trinuc_freq


def main(stayed_fasta, bg_fasta):
    seqs = {
        "locus_tag": [],
        "NT": [],
        "Content": [],

    }
    for seq_record in SeqIO.parse(stayed_fasta, "fasta"):
        res = trinucleotide_frequency(str(seq_record.seq))
        for trinuc, freq in res.items():
            seqs["locus_tag"].append(seq_record.id)
            seqs["NT"].append(trinuc)
            seqs["Content"].append(freq)
    seqs = pd.DataFrame(seqs)
    bg_seqs = {
        "locus_tag": [],
        "NT": [],
        "Content": [],

    }
    for seq_record in SeqIO.parse(bg_fasta, "fasta"):
        res = trinucleotide_frequency(str(seq_record.seq))
        for trinuc, freq in res.items():
            bg_seqs["locus_tag"].append(seq_record.id)
            bg_seqs["NT"].append(trinuc)
            bg_seqs["Content"].append(freq)


    bg_seqs = pd.DataFrame(bg_seqs)
    seqs["Background"] = False
    bg_seqs["Background"] = True

    data = pd.concat((seqs, bg_seqs))
    d = {
        "triNT": [],
        "pval": [],
    }
    for trinuc in data["NT"].unique():
        t = data[(data.NT == trinuc) & data.Background == False]["Content"]
        bg = data[(data.NT == trinuc) & data.Background == True]["Content"]
        stat, pval = mannwhitneyu(t, bg)
        d["triNT"].append(trinuc)
        d["pval"].append(pval)
    d = pd.DataFrame(d)
    _, pvals, _, _ = multipletests(d["pval"], method="fdr_bh")
    d["p.adjust"] = pvals
    print(d[d["p.adjust"] <= 0.05])
    data = data[data["NT"].isin(d[d["p.adjust"] <= 0.05]["triNT"])]
    fig = px.violin(
        data,
        x="NT",
        y="Content",
        color="Background",
        #points="all",
        box=True
    )
    fig.show()

if __name__ == '__main__':
    main(
        snakemake.input.stayed,
        snakemake.input.bg
    )