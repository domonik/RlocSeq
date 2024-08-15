from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd



def extract_fasta(tsv, type, output):
    df = pd.read_csv(tsv, sep="\t")
    stayed = df[df[type] == True]
    stayed = stayed.apply(lambda row: SeqRecord(Seq(row["RNASequence"]), id=row["locus_tag"], description=row["locus_tag"]), axis=1).tolist()
    SeqIO.write(stayed, output, format="fasta")


if __name__ == '__main__':
    extract_fasta(
        snakemake.input.tsv,
        snakemake.params.type,
        snakemake.output.fasta
    )