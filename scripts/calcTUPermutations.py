import time

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from multiprocessing import Pool
from plotly.subplots import make_subplots
from scripts.plotly_template import register_template

register_template()


def calc_mean_fraction(df, val):
    grouped = df.groupby('shortest_TU')[val].agg(
        fraction=lambda x: x.sum() / len(x),
        count='size'
    ).reset_index()
    grouped = grouped[grouped['count'] > 1]
    original_value = grouped["fraction"].mean()
    return original_value


def mp_wrapper(df):
    perm = np.random.permutation(len(df))
    jdf = df[["minusPMembran", "minusPCytosol"]]

    shuffled = df[['shortest_TU', 'TU_size']].iloc[perm].reset_index(drop=True)
    combined = pd.concat([shuffled.reset_index(drop=True), jdf.reset_index(drop=True)], axis=1)

    value1 = calc_mean_fraction(combined, "minusPMembran")
    value2 = calc_mean_fraction(combined, "minusPCytosol")
    return value1, value2


def main(df, cfg, outfile, num_permutations=1000, num_processes = 1):
    df = df[df["mRNA"]]
    original_value_mem = calc_mean_fraction(df, 'minusPMembran')
    original_value_cyt = calc_mean_fraction(df, 'minusPCytosol')

    df = df[["minusPMembran", "minusPCytosol", 'shortest_TU', 'TU_size']].copy()

    s = time.time()
    with Pool(num_processes) as p:
        values = p.starmap(mp_wrapper, [(df, ) for _ in range(num_permutations)])
    e = time.time()
    values = np.array(values)
    np.save(outfile, values)



if __name__ == '__main__':
    if "snakemake" in globals():
        tsv = snakemake.input.tsv
        cfg = snakemake.config
        num_permutations = snakemake.params.num_permutations
        threads = snakemake.threads
        outfile = snakemake.output.values
    else:
        tsv = "Pipeline/Paper/Supplementary/SummaryTable.tsv"
        import yaml
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        threads = 4
        num_permutations = 100
    df = pd.read_csv(tsv, sep="\t")
    main(df, outfile=outfile, num_processes=threads, cfg=cfg, num_permutations=num_permutations)