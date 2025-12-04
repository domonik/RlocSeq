import pandas as pd



def main(untreated, puro, tc, stayed, output):
    untreated = pd.read_csv(untreated, sep="\t")
    puro = pd.read_csv(puro, sep="\t")
    tc = pd.read_csv(tc, sep="\t")
    stayed = pd.read_csv(stayed, sep="\t")
    untreated["Dataset"] = "Untreated"
    puro["Dataset"] = "Puromycin"
    tc["Dataset"] = "Total cell"
    stayed["Dataset"] = "Untreated & Puromycin at Membrane"
    map_ = {
        "Untreated": "Membrane up",
        "Puromycin": "Membrane up",
        "Total cell": "Puromycin up",
        "Untreated & Puromycin at Membrane": "Membrane up",
    }
    df = pd.concat([untreated, puro, tc, stayed], ignore_index=True)
    df["testes subset"] = df["Dataset"].map(map_)
    other_cols = ["Dataset", "testes subset"]
    df = df[other_cols + [col for col in df.columns if col not in other_cols]]
    df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input.untreated,
        snakemake.input.puro,
        snakemake.input.tc,
        snakemake.input.stayed,
        snakemake.output.tsv,
    )