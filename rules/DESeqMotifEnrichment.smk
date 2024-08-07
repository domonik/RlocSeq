import os


rule getFasta:
    input:
        gff = config["gff"],
        fasta = config["fasta"]
    conda:
        "../envs/bedtools.yml"
    output:
        fasta = temporary(os.path.join(config["RUN_DIR"],"fasta/sequences.fa")),
    shell:
        "bedtools getfasta -s -fi {snakemake_input.fasta} -bed {snakemake_input.gff} > {output.fasta}"


rule extractFasta:
    input:
        deseq_result = config["DESeqResultTable"],
        gff = config["gff"],
        fasta = rules.getFasta.output.fasta
    conda:
        "../envs/biopython.yml"
    output:
        up_fasta = os.path.join(config["RUN_DIR"],"fasta/up.fa"),
        bg_fasta = os.path.join(config["RUN_DIR"],"fasta/bg.fa")
    script: "../scripts/extractFasta.py"

rule memeMotif:
    input:
        up_fasta = rules.extractFasta.output.up_fasta,
        bg_fasta = rules.extractFasta.output.bg_fasta,
    conda:
        "../envs/meme.yml"
    output:
        dir=directory(os.path.join(config["RUN_DIR"],"MEME/Motif_up/"))
    shell:
        """
        mkdir -p {output.dir}
        streme --p {input.up_fasta} --o {output.dir} --n {input.bg_fasta} --evalue --minw 5 --align left || true
        """

rule allMotif:
    default_target: True,
    input:
        file = rules.memeMotif.output