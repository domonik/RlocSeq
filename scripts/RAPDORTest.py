import pandas as pd
import time
from RAPDOR.datastructures import RAPDORData

s = time.time()

def cb(perc):
    global s
    e = time.time()

    print(perc, e-s)
    s= e


if __name__ == '__main__':
    design = pd.read_csv("../rapdor_test_design.tsv", sep="\t")
    df = pd.read_csv("../rapdor_test_intensities.tsv", sep="\t")
    data = RAPDORData(df=df, design=design, control="random_untreated")
    data.run_preprocessing()
    pvals, distribution = data.calc_anosim_p_value(100000, threads=1, callback=cb)
    data.determine_peaks(beta=10000)
    data.to_json("rapdortestfile.json")
    i = 0