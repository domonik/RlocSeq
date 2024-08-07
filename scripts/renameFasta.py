from Bio import SeqIO
from Bio.Seq import Seq
from DEplots import gff_helper
import pandas as pd




def main(stayed_gff, stayed_fasta, bg_gff, bg_fasta, stayed_out, bg_out, all_out):
    gff = gff_helper.read_gff3(stayed_gff)
    recs = []
    all_seqs = []
    for idx, seq_record in enumerate(SeqIO.parse(stayed_fasta, "fasta")):
        loc = gff.iloc[idx]
        seq_record.description = loc.locus_tag
        seq_record.id = loc.locus_tag
        seq_record.seq = Seq(str(seq_record.seq).replace("T", "U"))
        recs.append(seq_record)
        all_seqs.append(seq_record)
    SeqIO.write(recs, stayed_out, "fasta")
    gff = gff_helper.read_gff3(bg_gff)
    recs = []
    for idx, seq_record in enumerate(SeqIO.parse(bg_fasta, "fasta")):
        loc = gff.iloc[idx]
        seq_record.description = loc.locus_tag
        seq_record.id = loc.locus_tag
        seq_record.seq = Seq(str(seq_record.seq).replace("T", "U"))

        recs.append(seq_record)
        all_seqs.append(seq_record)
    SeqIO.write(recs, bg_out, "fasta")
    SeqIO.write(all_seqs, all_out, "fasta")






if __name__ == '__main__':
    main(
        snakemake.input.stayed_gff,
        snakemake.input.stayed,
        snakemake.input.bg_gff,
        snakemake.input.bg,
        snakemake.output.stayed,
        snakemake.output.bg,
        snakemake.output.all
    )