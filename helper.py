from multiprocessing import Pool
import multiprocessing
import random
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

multiprocessing.set_start_method("spawn", force=True)


def helper(df, threads):
    calls = [[df, x] for x in range(10000)]
    if threads > 1:
        with Pool(threads) as pool:
            result = pool.starmap(mp_wrapper, calls)
    else:
        result = [mp_wrapper(*call) for call in calls]
    return result


def mp_wrapper(df, x):
    padjs = flat_list = df['padjs'].sum()
    lfcs = df["lFCs"].sum()

    random.seed(x)
    random.shuffle(padjs)
    random.seed(x)
    random.shuffle(lfcs)
    nrGenes = df["numGenes"].to_numpy()
    df["cumsumGenes"] = np.cumsum(nrGenes)
    df["padjs"] = df[["cumsumGenes", "numGenes"]].apply(lambda row: padjs[row["cumsumGenes"] - row["numGenes"]: row[
        "cumsumGenes"]],axis=1)
    df["lFCs"] = df[["cumsumGenes", "numGenes"]].apply(lambda row: lfcs[row["cumsumGenes"] - row["numGenes"]: row[
        "cumsumGenes"]],axis=1)

    df["anySig+"] = df[["padjs", "lFCs"]].apply(lambda datas: any([pval <= .05 and lfc > 0 for pval, lfc in
                                                                   zip(datas["padjs"],datas["lFCs"])]),axis=1)
    df["anySig-"] = df[["padjs", "lFCs"]].apply(lambda datas: any([padj <= .05 and lfc < 0 for padj, lfc in
                                                                   zip(datas["padjs"],datas["lFCs"])]),axis=1)
    df["numSig+"] = df[["padjs", "lFCs"]].apply(lambda datas: sum([padj <= .05 and lfc > 0 for padj, lfc in
                                                                   zip(datas["padjs"],datas["lFCs"])]),axis=1)
    df["numSig-"] = df[["padjs", "lFCs"]].apply(lambda datas: sum([padj <= .05 and lfc < 0 for padj, lfc in
                                                                   zip(datas["padjs"],datas["lFCs"])]),axis=1)
    df["num+"] = df[["padjs", "lFCs"]].apply(lambda datas: sum([lfc > 0 for lfc in datas["lFCs"]]),axis=1)
    df["num-"] = df[["padjs", "lFCs"]].apply(lambda datas: sum([lfc < 0 for lfc in datas["lFCs"]]),axis=1)
    df["numGenes"] = df["lFCs"].apply(lambda x: len(x))

    m = df["numGenes"].max()
    d_up = {
        "TU-Size": []
    }
    d_down = {
        "TU-Size": []

    }
    for x in range(0,m + 1):
        d_up[f"{x} sig"] = []
        d_down[f"{x} sig"] = []
    for tu_size in range(1,m + 1):
        sub_df = df[df["numGenes"] == tu_size]
        d_up["TU-Size"].append(tu_size)
        d_down["TU-Size"].append(tu_size)
        for n_sig in range(0,m + 1):
            if n_sig > tu_size:
                n = np.nan
                n_down = np.nan
            else:
                n = sub_df[sub_df["numSig+"] == n_sig].shape[0]
                n_down = sub_df[sub_df["numSig-"] == n_sig].shape[0]
            d_up[f"{n_sig} sig"].append(n)
            d_down[f"{n_sig} sig"].append(n_down)
    out_df_up = pd.DataFrame(d_up)
    out_df_down = pd.DataFrame(d_down)
    m_up = out_df_up.iloc[0:, 1:].to_numpy()
    indices = np.arange(m_up.shape[0])
    tri_up = np.cumsum(m_up[indices, indices + 1][::-1])[::-1]
    whole_up = np.cumsum(np.nansum(m_up,axis=1)[::-1])[::-1] - out_df_up["0 sig"]
    out_df_up.loc[:, "greaterEqualTUSizeAllSignificant"] = tri_up
    out_df_up.loc[:, "greaterEqualTUSizeNumberOperons"] = whole_up
    out_df_up.loc[:, "greaterEqualTUSizePercentAllSignificant"] = np.round((tri_up / whole_up) * 100,decimals=1)
    return out_df_up.iloc[1]["greaterEqualTUSizePercentAllSignificant"]

def compare_results(files, sep, condition, base, padj_cutoff, lfc_cutoff, colors = ("rgb(0, 142, 151) ", "rgb(252, 76, 2)", "#2CA02C"), y: str = "Name"):
    fig = go.Figure()
    df = {"Enriched in": [], "Percentage": [], y: []}
    for name, file in files:
        df1 = pd.read_csv(file, sep=sep)
        df1 = df1[df1["baseMean"] > 0]
        df1_up = df1[(df1["padj"] <= padj_cutoff) & (df1["log2FoldChange"] >= lfc_cutoff)]
        df1_down = df1[(df1["padj"] <= padj_cutoff) & (df1["log2FoldChange"] <= -lfc_cutoff)]
        lup = len(df1_up)
        ldown = len(df1_down)
        lges = len(df1)

        df["Percentage"] += [(ldown/lges) * 100, ((lges - lup -ldown) / lges) * 100, (lup/lges) * 100]
        df[y] += [name, name, name]
        df["Enriched in"] += [base, "not Enriched", condition]
    df = pd.DataFrame(df)
    fig = px.bar(df, x="Percentage", y=y, color="Enriched in", color_discrete_sequence=colors, text="Percentage")
    fig.update_traces(marker_line_color='black',
                      marker_line_width=1)
    fig.update_traces(texttemplate='%{text:.1f}')
    fig.update_xaxes(range=[0, 100])

    return fig


def generateRDPMSpecInput(file, sep, annotation):
    anno = pd.read_csv(annotation, sep=",", index_col=0)

    treatments = []
    fractions = []
    replicates = []
    names = []
    drop = []
    tables = []
    df1 = pd.read_csv(file, sep=sep)
    df1 = df1.drop(df1[df1.index.str.contains("UTR")].index)
    for col in df1.columns:
        treat, fraction, rep = col.split("_")
        if fraction == "TC":
            drop.append(col)
        else:
            fraction = "Cytosol" if fraction == "C" else "Membrane"
            treat = "Puromycin+" if treat == "plusP" else "Puromycin-"
            treatments.append(treat)
            fractions.append(fraction)
            replicates.append(rep)
            names.append(col)
    df = pd.merge(anno, df1, left_index=True, right_index=True, how="right")
    print(df[df["gene_name"].str.contains("psbA2") == True])
    df = df.drop(drop, axis=1)
    df = df.rename({"gene_name": "Gene"}, axis=1)
    df.index.name = "locus_tag"
    design = {"Name": names, "Treatment": treatments, "Fraction": fractions, "Replicate": replicates}
    design = pd.DataFrame(design)
    return design, df
