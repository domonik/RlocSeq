import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from scipy.stats import mannwhitneyu
from scripts.plotly_template import register_template

register_template()

def main(snakemake_input, config, output):
    df = pd.read_csv(snakemake_input.rbp3, sep="\t")
    df = df[ [col for col in df.columns if col != "Localization"]]
    _map = {col: f"Rbp3_{col}" for col in df.columns if col != "locus_tag"} | {"locus_tag": "locus_tag"}

    df.columns = [_map[col] for col in df.columns]
    df = df[["locus_tag"] + [col for col in df.columns if col != "locus_tag"]]
    joined = pd.read_csv(snakemake_input.joined, sep="\t")
    joined = joined.merge(df, on="locus_tag", how="left")
    joined["Rbp3_Target"] = (~pd.isna(joined["Rbp3_log2FC"])) == True
    joined["Rbp3_Target"] = joined["Rbp3_Target"].astype("boolean")
    joined.loc[~(joined["mRNA"] == True), "Rbp3_Target"] = pd.NA
    joined.to_csv(output.useful_tsv, sep="\t", index=False)
    joined_out = joined[joined["mRNA"]]
    loc = ["Cellular thylakoid membrane"]
    cols = []
    for l in loc:
        col = f"Rbp3 Target and {l}".replace(" ", "_")
        cols.append(col)

        joined_out[col] = (joined_out["Rbp3_Target"] == True) & (joined_out[l] == True)
    moves = ["minusPMembran", "plusPMembran"]
    rna_keys = [key for key in list(config["RNA_TYPE_COLORS"].keys()) if key in joined_out.columns]
    subset = ["locus_tag", "log2FoldChange_plusPuromycin", "log2FoldChange_minusPuromycin", "Rbp3_Target"] + cols + moves + loc + rna_keys
    joined_out = joined_out[subset]
    joined_out.to_csv(output.trash_tsv_for_R, sep="\t", index=False)

    #fig.show()


if __name__ == '__main__':
    main(
        snakemake.input,
        snakemake.config,
        snakemake.output
    )
