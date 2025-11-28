import pandas as pd
from scipy.stats import ttest_rel, ttest_ind
from statsmodels.stats.multitest import multipletests


def calc_summary_and_pval(spots_df, probes):
    spots_df = spots_df[spots_df.Treatment != "Kasugamycin"]
    spots_df["subtractedIntensity"] = spots_df["ProbeIntensity"] - spots_df["background"]

    # Filter dataframe for selected probe
    df = spots_df[spots_df["Probe"].isin(probes)]
    df = df.groupby(["Treatment", "Probe", "Figure", ])["subtractedIntensity"].median().reset_index()
    # Run Welch’s t-test for each probe
    results = []
    for probe in df['Probe'].unique():
        untreated = df.query("Treatment == 'untreated' and Probe == @probe")['subtractedIntensity']
        puromycin = df.query("Treatment == 'Puromycin' and Probe == @probe")['subtractedIntensity']
        if len(untreated) == len(puromycin):
            tstat, pval = ttest_rel(untreated, puromycin, )
        else:
            tstat, pval = ttest_ind(untreated, puromycin, equal_var=False)
        results.append({'Probe': probe, 't-statistic': tstat, 'p-value': pval, "untreated": untreated.tolist(), "Puromycin": puromycin.tolist(), "untreated_mean": untreated.mean(), "puromycin_mean": puromycin.mean()})

    results_df = pd.DataFrame(results)
    reject, pvals_corr, _, _ = multipletests(results_df['p-value'], method='fdr_bh')

    results_df['adj p-value (FDR)'] = pvals_corr
    results_df['significant (FDR<0.05)'] = reject
    return df, results_df

def main(tsv, intensities_tsv, dist_out, intensities_out):

    df = pd.read_csv(tsv, sep="\t")
    df = df.rename({"Figure": "Replicate"}, axis=1)

    id_vars = ["Mask", "Image", "Probe", "Replicate", "Treatment"]
    mode = "mean"
    cols = [col for col in df.columns if mode in col and "bg_mean" not in col]
    d = {col: col.split("_")[1] for col in cols}
    df = df.rename(d, axis=1)
    cols = [p for p in d.values() if int(p) <= 20]
    df = df[cols + id_vars]

    df[cols] = (df[cols].values -  df[cols].values.min(axis=1)[:, None])/ (df[cols].values.max(axis=1)[:, None] -  df[cols].values.min(axis=1)[:, None])

    out_df = df.groupby(["Probe", "Replicate", "Treatment"])[[str(i) for i in range(20)]].mean()
    out_df["Peak position"] = out_df.idxmax(axis=1)

    out_df["n_spots"] = df.groupby(["Probe", "Replicate", "Treatment"])[str(0)].count()

    intensities = pd.read_csv(intensities_tsv, sep="\t")
    intensities, res_df = calc_summary_and_pval(spots_df=intensities, probes=["psbA2", "cupA", "fusB", "psaA", "16S rRNA"])
    out_df = out_df.reset_index()
    out_df = out_df[out_df["Treatment"] != "Kasugamycin"]
    map_ = {str(i): f"{i*26} nm" for i in range(20)}
    out_df = out_df.rename(map_, axis=1)
    out_df.to_csv(dist_out, sep="\t", index=False)
    res_df.to_csv(intensities_out, sep="\t", index=False)

if __name__ == '__main__':
    tsv = snakemake.input.tsv
    spots_df = snakemake.input.intensities

    main(tsv, spots_df, snakemake.output.dist, snakemake.output.intensities)