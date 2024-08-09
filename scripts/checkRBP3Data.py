import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from scipy.stats import mannwhitneyu


def main(snakemake_input, config, output):
    df = pd.read_csv(snakemake_input.rbp3, sep="\t")
    df = df[~df.Gene.str.contains("’")]
    joined = pd.read_csv(snakemake_input.joined, sep="\t")
    joined = joined.merge(df, on="locus_tag", how="left")
    joined["RBP3_Target"] = (~pd.isna(joined["log2FC"])) == True
    joined_out = joined[~(joined["UTR"] | joined["antisense"])]  # filter because they are not in the rbp3 dataset
    p=0
    loc = ["Cytoplasm", "Cell membrane", "Cellular thylakoid membrane"]
    cols = []
    for l in loc:
        col = f"RBP3 Target and {l}".replace(" ", "_")
        cols.append(col)

        joined_out[col] = (joined_out["RBP3_Target"] == True) & (joined_out[l] == True)
    moves = ["minusPMembran", "plusPMembran"]
    subset = ["locus_tag", "log2FoldChange_plusPuromycin", "log2FoldChange_minusPuromycin", "RBP3_Target"] + cols + moves + loc

    joined_out.to_csv(output.useful_tsv, sep="\t", index=False)
    joined_out = joined_out[subset]
    joined_out.to_csv(output.trash_tsv_for_R, sep="\t", index=False)
    data = joined[~pd.isna(joined["log2FoldChange_plusPuromycin"])]
    data = data[data["baseMean_plusPuromycin"] >= 100]
    pp = data[["RBP3_Target", "log2FoldChange_plusPuromycin"]]
    pp = pp.rename({"log2FoldChange_plusPuromycin": "log2FC"}, axis=1)
    pp = pp.rename({"log2FoldChange_plusPuromycin": "log2FC"}, axis=1)
    pp["Puromycin"] = "plus"
    mp = data[["RBP3_Target", "log2FoldChange_minusPuromycin"]]
    mp = mp.rename({"log2FoldChange_minusPuromycin": "log2FC"}, axis=1)
    mp["Puromycin"] = "minus"

    data = pd.concat((mp, pp))
    fig = go.Figure()

    fig.add_trace(
        go.Violin(
            x=data[data["RBP3_Target"]]["Puromycin"],
            y=data[data["RBP3_Target"]]["log2FC"],
            points="outliers",
            pointpos=0,
            box_visible=True,
            meanline_visible=True,
            legendgroup='Yes', scalegroup='Yes', name='Yes',
            side="positive"

        )
    )
    fig.add_trace(
        go.Violin(
            x=data[~data["RBP3_Target"]]["Puromycin"],
            y=data[~data["RBP3_Target"]]["log2FC"],
            points="outliers",
            pointpos=0,
            box_visible=True,
            meanline_visible=True,
            legendgroup='No', scalegroup='No', name='No',
            side="negative"

        )
    )
    fig.update_yaxes(
        title=dict(text="log<sub>2</sub>Fold Change"),

    )
    fig.update_xaxes(title=dict(text="Puromycin"))
    fig.update_traces(meanline_visible=True, width=0.5)
    fig.update_layout(violingap=0, violinmode='overlay', legend=dict(title=dict(text="RBP3_target")))
    res = mannwhitneyu(data[data["RBP3_Target"] & (data["Puromycin"] == "plus")]["log2FC"],
                       data[~data["RBP3_Target"] & (data["Puromycin"] == "plus")]["log2FC"])
    res2 = mannwhitneyu(data[data["RBP3_Target"] & (data["Puromycin"] == "minus")]["log2FC"],
                        data[~data["RBP3_Target"] & (data["Puromycin"] == "minus")]["log2FC"])
    fig.add_annotation(
        text=f"p-value: {res.pvalue:.2e}",
        x="plus",
        y=3.25,
        showarrow=False
    )
    fig.add_annotation(
        text=f"p-value: {res2.pvalue:.2e}",
        x="minus",
        y=3.25,
        showarrow=False

    )
    fig.update_xaxes(type='category')
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(
        config["RBP3A"]["layout"]
    )

    joined = joined[~pd.isna(joined["log2FC"])]
    df = joined
    fig = go.Figure()
    for i in ("Cytoplasm", "Membrane", "Cell membrane", "Cell inner membrane", "Cellular thylakoid membrane"):
        sdf = df[df[i] == True]
        ssdf = sdf[sdf["stayedAtMembrane"] == True]
        nsdf = sdf[~(sdf["stayedAtMembrane"] == True)]
        fig.add_trace(
            go.Scatter(
                x=ssdf["log2FC"],
                y=ssdf["log2FoldChange_plusPuromycin"],
                mode="markers",
                text=ssdf["gene_name"],
                marker=dict(line=dict(color="black", width=1)),
                name=i + " stayed"
            )
        )
        fig.add_trace(
            go.Scatter(
                x=nsdf["log2FC"],
                y=nsdf["log2FoldChange_plusPuromycin"],
                mode="markers",
                text=nsdf["gene_name"],
                name=i,
            )
        )
    sdf = df[pd.isna(df["Cytoplasm"])]
    fig.add_trace(
        go.Scatter(
            x=sdf["log2FC"],
            y=sdf["log2FoldChange_plusPuromycin"],
            mode="markers",
            text=sdf["gene_name"],
            name="Not annotated"
        )

    )
    fig.write_html(output.html)
    #fig.show()


if __name__ == '__main__':
    main(
        snakemake.input,
        snakemake.config,
        snakemake.output
    )