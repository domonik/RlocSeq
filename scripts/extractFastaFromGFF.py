from Bio import SeqIO
from DEplots import gff_helper
import pandas as pd



def extract_fasta(gff_file, tsv, type):
    df = pd.read_csv(tsv, sep="\t")
    df = df[df["mRNA"] == True]
    print(df)
    #stayed = df[(df["plusPMembran"] == True) & (df["Cellular thylakoid membrane"])]
    stayed = df[(df["stayedAtMembrane"] == True) & (df[type] == True)]
    bg = df[~((df["stayedAtMembrane"] == True))]
    #bg = df[~((df["plusPMembran"] == True))]


    gff = gff_helper.read_gff3(gff_file)
    stayed_gff = gff[gff["locus_tag"].isin(stayed["locus_tag"])]
    bg_gff = gff[gff["locus_tag"].isin(bg["locus_tag"])]
    stayed_gff["locus_tag"] = "locus_tag=" + stayed_gff["locus_tag"]
    bg_gff["locus_tag"] = "locus_tag=" + bg_gff["locus_tag"]

    return stayed_gff, bg_gff


if __name__ == '__main__':
    s, bg = extract_fasta(
        snakemake.input.gff,
        snakemake.input.tsv[0],
        type=snakemake.wildcards.type.replace("_", " ")
    )
    s.to_csv(snakemake.output.stayed, sep="\t", header=False, index=False)
    bg.to_csv(snakemake.output.bg, sep="\t", header=False, index=False)
