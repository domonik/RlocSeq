from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from DEplots import gff_helper
import pandas as pd



def extract_fasta(tsv, type, stayed_output, bg_output):
    df = pd.read_csv(tsv, sep="\t")


    #stayed = df[(df["plusPMembran"] == True) & (df["Cellular thylakoid membrane"])]
    stayed = df[(df["stayedAtMembrane"] == True) & (df[type])]
    stayed = stayed.apply(lambda row: SeqRecord(Seq(row["RNASequence"]), id=row["locus_tag"]), axis=1).tolist()
    SeqIO.write(stayed, stayed_output, format="fasta")
    bg = df[~((df["stayedAtMembrane"] == True))]
    bg = bg.apply(lambda row: SeqRecord(Seq(row["RNASequence"]), id=row["locus_tag"]), axis=1).tolist()
    SeqIO.write(bg, bg_output, format="fasta")


if __name__ == '__main__':
    extract_fasta(
        snakemake.input.tsv,
        snakemake.wildcards.type.replace("_", " "),
        snakemake.output.stayed,
        snakemake.output.bg,
    )

