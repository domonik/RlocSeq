import pandas as pd


def main(tsv, output):
    df = pd.read_csv(tsv, sep="\t")
    columns = [
        "locus_tag",
        "gene_name",
        "mRNA",
        "UTR",
        "asRNA",
        "ncRNA",
        "Cytoplasm",
        "Membrane",
        "Plasmamembrane",
        "Cellular thylakoid membrane",
        "TU",
        "Rbp3_Target"
    ]
    df =df[columns]
    df = df.rename({"gene_name": "Gene", "Cellular thylakoid membrane": "Thylakoid membrane"}, axis=1)
    df = df.set_index("locus_tag")
    df.to_csv(output, sep= "\t")

if __name__ == '__main__':
    if "snakemake" in globals():
        tsv = snakemake.input.tsv
        output = snakemake.output.tsv
    else:
        tsv = "Pipeline/Paper/Supplementary/SummaryTable.tsv"
        output = "foo.tsv"

    main(tsv, output)