import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def main(fasta, tsv_file, output_file):
    df = pd.read_csv(tsv_file, sep="\t")
    seqs = {
        "locus_tag": [],
        "RNASequence": []
    }
    for seq_record in SeqIO.parse(fasta, "fasta"):
        seqs["locus_tag"].append(seq_record.description)
        seqs["RNASequence"].append(str(seq_record.seq))
    seqs = pd.DataFrame(seqs)
    df = df.merge(seqs, on="locus_tag")
    df.to_csv(output_file, sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input.fasta,
        snakemake.input.tsv,
        snakemake.output.tsv
    )
