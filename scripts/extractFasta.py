import pandas as pd
from Bio import SeqIO



def extractFasta(deseq_result_file: str, gff_file, fasta_file, up_file, bg_file, config):
    df = pd.read_csv(deseq_result_file, sep="\t")

    df = df[df.index.str.contains("UTR")]
    up = df[df["log2FoldChange"] >= config["lfc_cutoff"]]
    up = up[up["padj"] <= config["padj_cutoff"]]
    gff = pd.read_csv(gff_file, sep="\t",
                      names=["chr", "feature", "feat", "start", "end", "dot", "strand", "score", "add"], skiprows=1)
    gff[config["add_column"]] = gff["add"].str.split(f"{config['add_column']}=").str[-1]
    gff["up"] = gff[config["add_column"]].isin(up.index)
    gff["match"] = gff["chr"] + ":" + (gff["start"] - 1).astype(str) + "-" + gff["end"].astype(str) + "(" + gff[
        "strand"] + ")"
    up = gff[gff["up"]]
    bg = gff[~gff["up"]]
    ups = []
    bgs = []
    d = up["match"].tolist()
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if str(seq_record.description) in d:
            ups.append(seq_record)
        else:
            bgs.append(seq_record)
    with open(up_file, "w") as handle:
        SeqIO.write(ups, handle, "fasta")

    with open(bg_file, "w") as handle:
        SeqIO.write(bgs, handle, "fasta")


if __name__ == '__main__':
    extractFasta(
        snakemake.input.deseq_result,
        snakemake.input.gff,
        snakemake.input.fasta,
        snakemake.output.up_fasta,
        snakemake.output.bg_fasta,
        snakemake.config
    )