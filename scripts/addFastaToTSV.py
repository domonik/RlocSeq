from Bio import SeqIO
from Bio.Seq import Seq
from DEplots import gff_helper
import pandas as pd




def main(gff, fasta, tsv, output_file):
    gff = gff_helper.read_gff3(gff)
    df = pd.read_csv(tsv, sep="\t")
    recs = []

    seqs = {
        "locus_tag": [],
        "RNASequence": []
    }
    for idx, seq_record in enumerate(SeqIO.parse(fasta, "fasta")):
        loc = gff.iloc[idx]
        assert loc["end"] - loc["start"] + 1 == len(seq_record.seq)
        seq_record.description = loc.locus_tag
        seq_record.id = loc.locus_tag
        seq_record.seq = Seq(str(seq_record.seq).replace("T", "U"))
        recs.append(seq_record)
        seqs["locus_tag"].append(seq_record.description)
        seqs["RNASequence"].append(str(seq_record.seq))
    seqs = pd.DataFrame(seqs)
    df = df.merge(seqs, on="locus_tag")
    df.to_csv(output_file, sep="\t", index=False)







if __name__ == '__main__':
    main(
        snakemake.input.gff,
        snakemake.input.fasta,
        snakemake.input.tsv,
        snakemake.output.tsv,
    )