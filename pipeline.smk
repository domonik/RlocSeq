import yaml
import pandas as pd
import numpy as np
import ast
import random
import os

NAMEMapping = {
    "C": "Cytosol",
    "M": "Membrane",
    "TC": "Total cell",
    "plusP": "+Puromycin",
    "minusP": "-Puromycin"
}



commit = "d976e73568037f58eaa16bd824f587781854cf36"

snakefile = github("domonik/SynDESeq2", path="snakefile.smk", commit=commit)
snakefile = "../SynDESeq2/snakefile.smk"
configfile: "configFile.yml"

config["RUNS"] = {}
with open(config["primary-file"], "r") as handle:
    config["RUNS"]["primary"] = yaml.safe_load(handle)
    config["primary"] = config["RUNS"]["primary"]

with open(config["secondary-file"], "r") as handle:
    config["RUNS"]["secondary"] = yaml.safe_load(handle)
    config["secondary"] = config["RUNS"]["secondary"]

with open(config["tertiary-file"], "r") as handle:
    config["sizeFactorEstimate"] = yaml.safe_load(handle)

with open(config["total_cell-file"], "r") as handle:
    config["total_cell"] = yaml.safe_load(handle)

with open(config["puromycinTC"], "r") as handle:
    config["puromycinTC"] = yaml.safe_load(handle)

with open(config["lightDarkLightConfig"], "r") as handle:
    config["RUNS"]["lightDarkLight"] = yaml.safe_load(handle)

with open(config["lightDarkDarkConfig"], "r") as handle:
    config["RUNS"]["lightDarkDark"] = yaml.safe_load(handle)

with open(config["lightDarkLightPuroConfig"], "r") as handle:
    config["RUNS"]["lightDarkLightPuro"] = yaml.safe_load(handle)

with open(config["lightDarkDarkPuroConfig"], "r") as handle:
    config["RUNS"]["lightDarkDarkPuro"] = yaml.safe_load(handle)

with open(config["lightDarkLightTCConfig"], "r") as handle:
    config["RUNS"]["lightDarkLightTC"] = yaml.safe_load(handle)

with open(config["lightDarkLightvsLightTCConfig"], "r") as handle:
    config["RUNS"]["lightDarkLightvsLightTC"] = yaml.safe_load(handle)

with open(config["lightDarkLightvsLightPlusPuroTCConfig"], "r") as handle:
    config["RUNS"]["lightDarkLightvsLightPlusPuroTC"] = yaml.safe_load(handle)





rule joinLandLDL:
    input:
        light = config["RUNS"]["primary"]["CountFile"],
        light_dark = config["RUNS"]["lightDarkLightTC"]["CountFile"],
    output:
        file = "Pipeline/BothDatasetsJoined.tsv"
    run:
        import pandas as pd
        import re
        df1 = pd.read_csv(input.light, sep=",")
        df2 = pd.read_csv(input.light_dark, sep=",")
        def rename_col(col):
            # Step 1: Replace 'minusP' with 'minP'
            col = col.replace('minusP','minP')

            # Step 2: Insert 'HL' between region and replicate
            # Pattern: (prefix)_(region)_(replicate)
            match = re.match(r"(.+P)_(\w+)_(R\d+)",col)
            if match:
                prefix, region, replicate = match.groups()
                return f"{prefix}_{region}_HL_{replicate}"
            else:
                return col


        df1.columns = [rename_col(col) for col in df1.columns]
        df1.index = df1.iloc[:, 0]
        df2.index = df2.iloc[:, 0]
        df1 = df1[[col for col in df1.columns[1:]]]
        df2 = df2[[col for col in df2.columns[1:]]]
        df_all = df1.merge(df2,left_index=True,right_index=True)
        df_all.to_csv(output.file, sep=",")



rule coLocateTranscripts:
    input:
        tus = "initialData/Supplementary_Information1_6803MS_Final_added_annotation.xlsx"
    output:
        tus =  "Pipeline/TUAnalysis/TUFile.tsv"
    run:
        import pandas as pd
        df = pd.read_excel(input.tus, skiprows=1)
        sdf = df[[ "Sense tags" , "TU ID"]]
        asdf = df[["Antisense tags", "Antisense TUs"]]
        asdf["Antisense TUs"] = asdf["Antisense TUs"].str.split(", ")
        #asdf["Antisense tags"] = asdf["Antisense tags"].str.split(", ")
        asdf = asdf.explode("Antisense TUs")

        #sdf["Sense tags"] = sdf["Sense tags"].str.split(", ")
        #sdf = sdf.explode("Sense tags")
        sdf = sdf[~sdf["TU ID"].isin(asdf["Antisense TUs"].unique())]
        sdf = sdf.rename({"Sense tags": "locus_tags"}, axis=1)
        asdf = asdf.rename({"Antisense tags": "locus_tags", "Antisense TUs": "TU ID"}, axis=1)
        sdf.loc[:, "Antisense"] = False
        asdf.loc[:, "Antisense"] = True
        df = pd.concat((sdf, asdf))
        df = df[~df["locus_tags"].isna()]
        df.to_csv(output.tus, sep="\t", index=False)





module other_workflow:
    snakefile: snakefile
    config: config["RUNS"]["primary"]

use rule * from other_workflow as other_*

module second_workflow:
    snakefile: snakefile
    config: config["RUNS"]["secondary"]

use rule * from second_workflow as second_*

module third_workflow:
    snakefile: snakefile
    config: config["sizeFactorEstimate"]

use rule * from third_workflow as third_*

module tc_workflow:
    snakefile: snakefile
    config: config["total_cell"]

use rule * from tc_workflow as tc_*


module MemeMotif:
    snakefile: "rules/DESeqMotifEnrichment.smk"
    config: config["memeMotif"]

use rule * from MemeMotif as meme_*

module MemeMotif2:
    snakefile: "rules/DESeqMotifEnrichment.smk"
    config: config["memeMotif2"]

use rule * from MemeMotif2 as meme2_*


module MemeMotifLDL:
    snakefile: "rules/DESeqMotifEnrichment.smk"
    config: config["memeMotifLD_L"]

use rule * from MemeMotifLDL as meme_ldl_*

module MemeMotifLDD:
    snakefile: "rules/DESeqMotifEnrichment.smk"
    config: config["memeMotifLD_D"]

use rule * from MemeMotifLDD as meme_ldd_*

module MemeMotifLDDPuro:
    snakefile: "rules/DESeqMotifEnrichment.smk"
    config: config["memeMotifLD_DPuro"]

use rule * from MemeMotifLDDPuro as meme_lddpuro_*

module MemeMotifLDLPuro:
    snakefile: "rules/DESeqMotifEnrichment.smk"
    config: config["memeMotifLD_LPuro"]

use rule * from MemeMotifLDLPuro as meme_ldlpuro_*

def checkLocusTagspvals(other_df, col):
    def check(ltags):
        l = ltags.split(", ")
        d = []
        for ltag in l:
            if ltag in other_df.index:
                d.append(other_df.loc[ltag, col])
        return d

    return check


rule create_uniprot_annotation:
    input:
        uniprot = rules.other_downloadOrganismGOTerms.output.go_terms,
        df = config["RUNS"]["primary"]["CountFile"]
    output:
        tsv = "Pipeline/uniprotMapping.tsv"
    run:
        import pandas as pd
        uprot = pd.read_csv(input.uniprot, sep="\t")
        uprot = uprot.assign(**{"Gene Names (ordered locus)": uprot["Gene Names (ordered locus)"].str.split(r"\s*[;/]\s*")})
        uprot = uprot.explode("Gene Names (ordered locus)").reset_index(drop=True)
        df = pd.read_csv(input.df, sep=config["RUNS"]["primary"]["initial_sep"])
        mapping = df.merge(uprot, left_on="Geneid", right_on="Gene Names (ordered locus)")
        mapping = mapping[["Geneid", "Entry"]]
        mapping.to_csv(output.tsv, sep="\t", index=False)

def any_sig(data):
    pvals = data["padjs"]
    lfcs = data["lFCs"]
    return any([pval < 0.05 and lfc > 0 for pval, lfc in zip(pvals, lfcs)])



rule checkTUs:
    input:
        tu_file = rules.coLocateTranscripts.output.tus,
        deseqRes = rules.other_extractDESeqResult.output.result_table
    output:
        file = "Pipeline/TUAnalysis/tulist_c{condition}_vs_b{baseline}.tsv"
    run:
        df = pd.read_csv(input.tu_file, sep="\t")
        de_res = pd.read_csv(input.deseqRes, sep="\t")
        df["padjs"] = df["locus_tags"].apply(checkLocusTagspvals(de_res, "padj"))
        df["lFCs"] = df["locus_tags"].apply(checkLocusTagspvals(de_res, "log2FoldChange"))
        df["anySig+"] = df[["padjs", "lFCs"]].apply(lambda datas:  any([pval <= .05 and lfc > 0.8 for pval, lfc in zip(datas["padjs"], datas["lFCs"])]), axis=1)
        df["anySig-"] = df[["padjs", "lFCs"]].apply(lambda datas:  any([padj <= .05 and lfc < -0.8 for padj, lfc in zip(datas["padjs"], datas["lFCs"])]), axis=1)
        df["numSig+"] = df[["padjs", "lFCs"]].apply(lambda datas:  sum([padj <= .05 and lfc > 0.8 for padj, lfc in zip(datas["padjs"], datas["lFCs"])]), axis=1)
        df["numSig-"] = df[["padjs", "lFCs"]].apply(lambda datas:  sum([padj <= .05 and lfc < -0.8 for padj, lfc in zip(datas["padjs"], datas["lFCs"])]), axis=1)
        df["num+"] = df[["padjs", "lFCs"]].apply(lambda datas:  sum([lfc > 0 for lfc in  datas["lFCs"]]), axis=1)
        df["num-"] = df[["padjs", "lFCs"]].apply(lambda datas:  sum([lfc < 0 for lfc in datas["lFCs"]]), axis=1)
        df["numGenes"] = df["lFCs"].apply(lambda x: len(x))
        df = df[~df["padjs"].apply(lambda x: bool(np.any(np.isnan(x))))] # filter TUs without values e.g. low counts
        df.to_csv(output.file, sep="\t", index=False)

def parse_list_string(x):
        try:
            return ast.literal_eval(x)
        except ValueError:
            return x




rule createBackground:
    input:
        df = rules.checkTUs.output.file
    output:
        tu_background = "Pipeline/TUAnalysis/TUBackground_c{condition}_vs_b{baseline}_up.tsv",
    threads: 6
    run:
        import pandas as pd
        from helper import helper

        import numpy as np
        import random
        df = pd.read_csv(input.df, sep="\t", converters={'values': eval})
        df['padjs'] = df['padjs'].apply(parse_list_string)
        df['lFCs'] = df['lFCs'].apply(parse_list_string)
        result = helper(df, threads=threads)
        df = pd.DataFrame({"randomPercantages": result})
        df.to_csv(output.tu_background, sep="\t", index=False)



rule createTUTable2:
    input:
        tsv=rules.checkTUs.output.file
    output:
        tu_table = "Pipeline/TUAnalysis/TUSummary2_c{condition}_vs_b{baseline}_up.tsv",
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_csv(input.tsv,sep="\t")
        m = df["numGenes"].max()
        d_up = {
            "TU-Size": []
        }
        df = df[df["anySig+"]]
        for x in range(0,m + 1):
            d_up[f"{x} up"] = []
        for tu_size in range(1,m + 1):
            sub_df = df[df["numGenes"] == tu_size]
            d_up["TU-Size"].append(tu_size)
            for n_sig in range(0,m + 1):
                if n_sig > tu_size:
                    n = np.nan
                    n_down = np.nan
                else:
                    n = sub_df[sub_df["num+"] == n_sig].shape[0]
                d_up[f"{n_sig} up"].append(n)
        out_df_up = pd.DataFrame(d_up)
        m_up = out_df_up.iloc[0:, 1:].to_numpy()
        indices = np.arange(m_up.shape[0])
        tri_up = np.cumsum(m_up[indices, indices + 1][::-1])[::-1]
        whole_up = np.cumsum(np.nansum(m_up,axis=1)[::-1])[::-1] - out_df_up["0 up"]
        out_df_up.loc[:, "greaterEqualTUSizeAllSignificant"] = tri_up
        out_df_up.loc[:, "greaterEqualTUSizeNumberOperons"] = whole_up
        out_df_up.loc[:, "greaterEqualTUSizePercentAllSignificant"] = np.round((tri_up / whole_up) * 100,decimals=1)
        out_df_up.to_csv(output.tu_table, sep="\t", index=False)


rule download_subcell:
    output:
        tsv = "Pipeline/SubcellularLoc/uniprot_anno.tsv"
    script: "scripts/downloadCellularLocation.py"


rule createTUTable:
    input:
        tsv = rules.checkTUs.output.file
    output:
        tu_table = "Pipeline/TUAnalysis/TUSummary_c{condition}_vs_b{baseline}_up.tsv",
        tu_table_down = "Pipeline/TUAnalysis/TUSummary_c{condition}_vs_b{baseline}_down.tsv",
    run:
        import pandas as pd
        import numpy as np
        df = pd.read_csv(input.tsv, sep="\t")
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
        whole_up = np.cumsum(np.nansum(m_up, axis=1)[::-1])[::-1] - out_df_up["0 sig"]
        out_df_up.loc[:, "greaterEqualTUSizeAllSignificant"] = tri_up
        out_df_up.loc[:, "greaterEqualTUSizeNumberOperons"] = whole_up
        out_df_up.loc[:, "greaterEqualTUSizePercentAllSignificant"] = np.round((tri_up / whole_up) * 100, decimals=1)

        m_down = out_df_down.iloc[0:, 1:].to_numpy()
        indices = np.arange(m_down.shape[0])
        tri_down = np.cumsum(m_down[indices, indices + 1][::-1])[::-1]
        whole_down = np.cumsum(np.nansum(m_down,axis=1)[::-1])[::-1] - out_df_down["0 sig"]
        out_df_down.loc[:, "greaterEqualTUSizeAllSignificant"] = tri_down
        out_df_down.loc[:, "greaterEqualTUSizeNumberOperons"] = whole_down
        out_df_down.loc[:, "greaterEqualTUSizePercentAllSignificant"] = np.round((tri_down / whole_down) * 100,decimals=1)

        out_df_up.to_csv(output.tu_table, sep="\t", index=False)
        out_df_down.to_csv(output.tu_table_down, sep="\t", index=False)




rule plotHistogramm:
    input:
        bg = rules.createBackground.output.tu_background,
        values = rules.createTUTable.output.tu_table
    output:
        svg = "Pipeline/TUAnalysis/TU_Histogram_c{condition}_vs_b{baseline}_up.svg",
    run:
        import pandas as pd
        import plotly.graph_objs as go
        bg = pd.read_csv(input.bg, sep="\t")
        df = pd.read_csv(input.values, sep="\t")
        comparePerc = df.iloc[1]["greaterEqualTUSizePercentAllSignificant"]
        print(comparePerc)
        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=bg["randomPercantages"],
            name="random Permutations"
        ))
        fig.add_vline(
            x=comparePerc,
            showlegend=True,
            name="greaterEqualTUSizePercentAllSignificant"
        )
        fig.show()


rule summarizeResult:
    input:
        deseqResMinus = rules.other_extractDESeqResult.output.result_table,
        deseqResPlus = rules.second_extractDESeqResult.output.result_table,
        mapping = rules.other_downloadOrganismGOTerms.output.go_terms,
        gff = "initialData/20210217_syne_onlyUnique_withFeat.gff3"
    output:
        file = "Pipeline/RNAType/RNASummary_c{condition}_vs_b{baseline}.tsv"
    run:
        import numpy as np
        df_minus = pd.read_csv(input.deseqResMinus, sep="\t")
        df_plus = pd.read_csv(input.deseqResPlus, sep="\t")
        gff = pd.read_csv(input.gff, sep="\t", names=["chr", "feature", "type", "start", "end", "dot", 	"strand", "score", "add_info"], skiprows=1)

        mapping = pd.read_csv(input.mapping, sep="\t")
        mapping.index = mapping["Gene Names (ordered locus)"]
        mapping = mapping.drop("Gene Names (ordered locus)", axis=1)
        df_minus = df_minus.merge(mapping, left_index=True, right_index=True, how="left")
        df_plus = df_plus.merge(mapping, left_index=True, right_index=True, how="left")
        out_d = {
            "RNA Type": [],
            "Total": [],
            "minus Puromycin Cytosol": [],
            "minus Puromycin Membran": [],
            "plus Puromycin Cytosol": [],
            "plus Puromycin Membran": [],
        }
        filters = {"mRNA": (r'\b(?:sl[rl]|ss[rl])\w*(?![-_])(?<!-as)\b', "Gene Names"), "5'UTR": ("5UTR", "idx"), "3'UTR": ("3UTR", "idx"), "tRNA": ("6803t", "idx"), "asRNA": ("-as", "idx"), "Photosystem": ("psa|psb", "Gene Names"), "Ribosome": ("rps|rpl", "Gene Names"), "CRISPR": ("CRISPR", "idx"), "ncRNA": ("ncr", "idx"), }
        for name, (filter, column) in filters.items():

            out_d["RNA Type"].append(name)
            for idx, (s, df) in enumerate({"minus Puromycin": df_minus, "plus Puromycin": df_plus}.items()):
                if column == "idx":
                    sub_df = df[df.index.str.contains(filter)]
                else:
                    sub_df = df[df[column].str.contains(filter) == True]
                if idx == 0:
                    out_d["Total"].append(sub_df.shape[0])

                condition = sub_df[(sub_df["log2FoldChange"] >= config["primary"]["log2FCCutOff"]) & (sub_df["padj"] <= config["primary"]["pAdjCutOff"])].shape[0]
                baseline = sub_df[(sub_df["log2FoldChange"] <= -config["primary"]["log2FCCutOff"]) & (sub_df["padj"] <= config["primary"]["pAdjCutOff"])].shape[0]
                out_d[f"{s} Membran"].append(condition)
                out_d[f"{s} Cytosol"].append(baseline)

        out_df = pd.DataFrame(out_d)
        arr = out_df.iloc[:, 2:].to_numpy()
        percent = np.round(( arr / out_df["Total"].to_numpy().reshape(arr.shape[0], 1)) * 100, decimals=1)
        percent = np.char.add(percent.astype(str), "%)")
        percent = np.char.add(" (", percent)
        out_df.iloc[:, 2:] = np.char.add(out_df.iloc[:, 2:].values.astype(str), percent)
        out_df.to_csv(output.file, sep="\t", index=False)

rule joinTables:
    input:
        deseqResMinus = rules.other_extractDESeqResult.output.result_table,
        deseqResPlus = rules.second_extractDESeqResult.output.result_table,
        add_data = "initialData/locusTagtoGeneName.csv",
        tu_info = rules.coLocateTranscripts.output.tus
    output:
        tsv = "Pipeline/Summary/joinedDESeq_c{condition}_vs_b{baseline}.tsv",
        output_for_R_because_this_stupid_language_cannot_even_read_tsv_files = "Pipeline/Summary/joinedDESeqForstupidR_c{condition}_vs_b{baseline}.tsv",
        stayed = "Pipeline/Summary/stayedAtMembrane_c{condition}_vs_b{baseline}.tsv",
    run:
        df_minus = pd.read_csv(input.deseqResMinus, sep="\t")
        df_plus = pd.read_csv(input.deseqResPlus, sep="\t")
        df = df_minus.join(df_plus, rsuffix="_plusPuromycin", lsuffix="_minusPuromycin", how="outer")
        df = df[~df.index.str.contains("6803r")]

        df["minusPMembran"] = (df["log2FoldChange_minusPuromycin"] >= config["primary"]["log2FCCutOff"]) & (df["padj_minusPuromycin"] <= config["primary"]["pAdjCutOff"])
        df["minusPCytosol"] = (df["log2FoldChange_minusPuromycin"] <= -config["primary"]["log2FCCutOff"]) & (df["padj_minusPuromycin"] <= config["primary"]["pAdjCutOff"])
        df["plusPCytosol"] = (df["log2FoldChange_plusPuromycin"] <= -config["primary"]["log2FCCutOff"]) & (df["padj_plusPuromycin"] <= config["primary"]["pAdjCutOff"])
        df["plusPMembran"] = (df["log2FoldChange_plusPuromycin"] >= config["primary"]["log2FCCutOff"]) & (df["padj_plusPuromycin"] <= config["primary"]["pAdjCutOff"])
        df["shiftMtoC"] = (df["minusPMembran"] & df["plusPCytosol"])
        df["shiftCtoM"] = (df["plusPMembran"] & df["minusPCytosol"])
        df["stayedAtMembrane"] = df["plusPMembran"] & df["minusPMembran"]
        df["stayedInCytosol"] = df["plusPCytosol"] & df["minusPCytosol"]
        df["UTR"] = df.index.str.contains("UTR")
        df["5UTR"] = df.index.str.contains("5UTR")
        df["3UTR"] = df.index.str.contains("3UTR")
        df["asRNA"] = df.index.str.contains("-as")
        df["ncRNA"] = df.index.str.contains("ncr|ncl")
        df["tRNA"] = df.index.str.contains("6803t")
        df["CRISPR"] = df.index.str.contains("CRISPR")
        df["mRNA"] = (df.index.str.startswith("sll") | df.index.str.startswith("slr") | df.index.str.startswith("ssl") | df.index.str.startswith("ssr") | df.index.str.startswith("sml" )| df.index.str.startswith("smr") | df.index.str.startswith("sgl") | df.index.str.startswith("Norf")) & ~(df["UTR"] | df["asRNA"])
        add_data = pd.read_csv(input.add_data, sep=",")

        df = df.merge(add_data, left_index=True, right_on="locus_tag")
        df.loc[df["5UTR"], "gene_name"] = df.loc[df["5UTR"], "gene_name"] + " 5UTR"
        df.loc[df["3UTR"], "gene_name"] = df.loc[df["5UTR"], "gene_name"] + " 3UTR"
        df["Photosystem mRNA"] = df.gene_name.str.contains("psa|psb")
        df["Ribosome mRNA"] = df.gene_name.str.contains("rpl|rps")
        df = df[df.columns[-2:].to_list() + df.columns[:-2].to_list()]
        tus = pd.read_csv(input.tu_info, sep="\t")

        tus_expanded = (
            tus.assign(locus_tag=tus.loc[:, "locus_tags"].str.split(", "))
            .explode("locus_tag")
            .rename(columns={tus.columns[1]: "TU"})
        )
        # Add TU column to df
        df["locus_tag_norm"] = df["locus_tag"].str.replace(r'_(5UTR|3UTR)$', '', regex=True)
        tus_expanded["locus_tag"] = tus_expanded["locus_tag"].str.strip()
        tus_expanded = tus_expanded[["locus_tag", "TU"]]
        tus_grouped = (
            tus_expanded.groupby("locus_tag")["TU"]
            .apply(lambda x: ",".join(sorted(set([str(i) for i in x if ~pd.isna(i)]))))
            .reset_index()
        )
        df = df.merge(
            tus_grouped.rename(columns={"locus_tag": "locus_tag_norm"}),
            on="locus_tag_norm",
            how="left"
        )
        #df = df.merge( tus_expanded.rename(columns={"locus_tag": "locus_tag_norm"}),on="locus_tag_norm",how="left")

        df.to_csv(output.tsv, sep="\t", index=False)

        rtkeys = [key for key in list(config["RNA_TYPE_COLORS"].keys()) if key in df.columns]
        df_for_r = df[["locus_tag", "stayedAtMembrane", "stayedInCytosol", "mRNA", "log2FoldChange_plusPuromycin", "log2FoldChange_minusPuromycin"] + rtkeys]

        df = df[df["stayedAtMembrane"] == True]

        df.to_csv(output.stayed, sep="\t", index=False)
        df_for_r.to_csv(output.output_for_R_because_this_stupid_language_cannot_even_read_tsv_files, sep="\t", index=False)


rule Table1:
    input:
        tsv = expand(rules.joinTables.output.tsv, condition=["M"], baseline=["C"])[0],
    output:
        csv = "Pipeline/Paper/Table1/Table1.csv"
    run:
        df = pd.read_csv(input.tsv, sep="\t")
        names = ["mRNA", "5UTR", "3UTR", "tRNA", "asRNA", "Photosystem mRNA", "Ribosome mRNA", "CRISPR", "ncRNA"]
        summary = [
            {
                "RNA type": "All",
                "minus Puromycin Total":  df[~(df["baseMean_minusPuromycin"] == 0)].shape[0],
                "minus Puromycin Cytosol": df["minusPCytosol"].sum(),
                "minus Puromycin Membrane": df["minusPMembran"].sum(),
                "plus Puromycin Total": df[~(df["baseMean_plusPuromycin"] == 0)].shape[0],
                "plus Puromycin Membrane": df["plusPMembran"].sum(),
                "plus Puromycin Cytosol": df["plusPCytosol"].sum(),
                "both Total":
                    df[~(df["baseMean_minusPuromycin"] == 0) & ~(df["baseMean_plusPuromycin"] == 0)].shape[
                        0],
                "both in Cytosol": df["stayedInCytosol"].sum(),
                "both at Membrane": df["stayedAtMembrane"].sum(),
            }
        ]

        for col in names:
            subset = df[df[col] == True]  # rows where this RNA type is True
            summary.append({
                "RNA type": col,
                "minus Puromycin Total":  subset[~(subset["baseMean_minusPuromycin"] == 0)].shape[0],
                "minus Puromycin Cytosol": subset["minusPCytosol"].sum(),
                "minus Puromycin Membrane": subset["minusPMembran"].sum(),
                "plus Puromycin Total": subset[~(subset["baseMean_plusPuromycin"] == 0)].shape[0],
                "plus Puromycin Membrane": subset["plusPMembran"].sum(),
                "plus Puromycin Cytosol": subset["plusPCytosol"].sum(),
                "both Total": subset[~(subset["baseMean_minusPuromycin"] == 0) & ~(subset["baseMean_plusPuromycin"] == 0)].shape[0],
                "both in Cytosol": subset["stayedInCytosol"].sum(),
                "both at Membrane": subset["stayedAtMembrane"].sum(),

            })

        out_df = pd.DataFrame(summary)
        for puro in ["minus", "plus"]:
            loc = [f"{puro} Puromycin Cytosol", f"{puro} Puromycin Membrane"]
            arr = out_df.loc[:, loc].to_numpy()
            percent = np.round(( arr / out_df[f"{puro} Puromycin Total"].to_numpy().reshape(arr.shape[0], 1)) * 100, decimals=1)
            percent = np.char.add(percent.astype(str), "%)")
            percent = np.char.add(" (", percent)
            out_df.loc[:, loc] = out_df.loc[:, loc].astype(str)
            out_df.loc[:, loc] = np.char.add(out_df.loc[:, loc].values.astype(str), percent)
        loc = [f"both in Cytosol", f"both at Membrane"]
        arr = out_df.loc[:, loc].to_numpy()
        percent = np.round((arr / out_df[f"both Total"].to_numpy().reshape(
            arr.shape[0],1)) * 100,decimals=1)
        percent = np.char.add(percent.astype(str),"%)")
        percent = np.char.add(" (",percent)
        out_df.loc[:, loc] = out_df.loc[:, loc].astype(str)
        out_df.loc[:, loc] = np.char.add(out_df.loc[:, loc].values.astype(str),percent)
        out_df.to_csv(output.csv, sep=",", index=False)


rule GOTermEnrichmentStayed:
    input:
        df = rules.joinTables.output.output_for_R_because_this_stupid_language_cannot_even_read_tsv_files,
        annotation_db = rules.other_generateOrgDB.output.annotation_db,
    params:
        cfg = config["primary"]
    output:
        up = "Pipeline/Summary/GO/joinedDESeq_c{condition}_vs_b{baseline}_ud_up.tsv",
        down = "Pipeline/Summary/GO/joinedDESeq_c{condition}_vs_b{baseline}_ud_down.tsv"
    conda:
        "envs/clusterprofiler.yml"
    script: "scripts/GOTermStayedAtMembrane.R"


rule GSEAGO:
    input:
        enrich1 = rules.other_GSEAGO.output.enriched,
        gsdata1 = rules.other_GSEAGO.output.gsdata,
        enrich2= rules.other_GSEAGO.output.enriched,
        gsdata2=rules.other_GSEAGO.output.gsdata,
    output:
        html = "Pipeline/Paper/GSEGOplasmaMembrane.html"
    script: "scripts"


rule extractMembraneAndSolubleProteins:
    input:
        membrane = "initialData/MembraneProteins.xlsx",
        soluble = "initialData/SolubleProteins.xlsx"
    output:
        tsv = "Pipeline/ProteinLocalization/proteins.tsv"
    run:
        from scipy.stats import ttest_ind, false_discovery_control

        membrane_df = pd.read_excel(input.membrane, skiprows=2)
        cols = ["Accession", "Gene_symbol"] + [col for col in membrane_df.columns if "Ctrl1" in col]
        membrane_df =membrane_df[cols ]
        membrane_df = membrane_df[~membrane_df["Accession"].isna()]

        membrane_df["Mean"] = membrane_df[cols[2:]].mean(axis=1)
        soluble_df = pd.read_excel(input.soluble, skiprows=2)
        soluble_df = soluble_df[cols]
        soluble_df["Mean"] = soluble_df[cols[2:]].mean(axis=1)

        soluble_df = soluble_df[~soluble_df["Accession"].isna()]

        stats = []




        df = membrane_df.merge(soluble_df, left_on="Accession", right_on="Accession", suffixes=["_membrane", "_soluble"])
        df = df.drop("Gene_symbol_membrane", axis=1).rename({"Gene_symbol_soluble": "Gene_symbol"}, axis=1)
        df = df.fillna(0)
        membrane = df[[col + "_membrane" for col in cols[2:]]].values
        soluble = df[[col + "_soluble" for col in cols[2:]]].values
        for row1, row2 in zip(membrane, soluble):
            try:
                stat, pval = ttest_ind(a=row1, b=row2, equal_var=False)
            except ValueError:
                stat = np.nan
            stats.append(pval)
        df["pval"] = stats

        df = df[~df["pval"].isna()]
        df["padj"] = false_discovery_control(df["pval"])
        df["log2FoldChange"] = np.log2(df["Mean_membrane"]) - np.log2(df["Mean_soluble"])
        df = df[~df['log2FoldChange'].isin([np.inf, -np.inf])]
        df.to_csv(output.tsv, sep="\t", index=False)


hovertemplate = 'Gene: <b>%{text}</b><br>RNA padj: %{customdata[0]:.2e}<br>Protein padj: %{customdata[1]:.2e}'


rule correlationRNAvsProteins:
    input:
        proteins = rules.extractMembraneAndSolubleProteins.output.tsv,
        rna = rules.second_extractDESeqResult.output.result_table,
        summary =  rules.joinTables.output.tsv
    output:
        html = "Pipeline/ProteinLocalization/correlation_c{condition}_vs_b{baseline}.html"
    run:
        import plotly.graph_objs as go
        from scipy.stats import spearmanr
        res = pd.read_csv(input.summary,sep="\t", index_col=0)
        res_ids = res[res["stayedAtMembrane"]].index

        rna_df = pd.read_csv(input.rna, sep="\t")
        rna_df["Mean_membrane"] = rna_df["baseMean"] * rna_df["log2FoldChange"]
        rna_df = rna_df.rename({"baseMean": "Mean_soluble"}, axis=1)

        protein_df = pd.read_csv(input.proteins, sep="\t")
        print(protein_df)
        df = rna_df.merge(protein_df, left_index=True, right_on = "Accession", suffixes=["_RNA", "_Protein"])
        df = df[~df["log2FoldChange_RNA"].isna() & ~df["log2FoldChange_Protein"].isna()]
        #df  = df[(df["padj_RNA"] <= 0.05) & (df["padj_Protein"] <= 0.05)]
        df["hover"] = df["Accession"] + "-" + df["Gene_symbol"]

        r, pval = spearmanr(df["log2FoldChange_RNA"], df["log2FoldChange_Protein"])
        print(r, pval)
        fig = go.Figure()
        cd = np.stack((df["padj_RNA"], df["padj_Protein"]),axis=1)
        fig.add_annotation(
            x=0,
            y=1.1,
            xref="x domain",
            yref="y domain",
            text=f"Spearman:\n\tR: {r:.2f}\n\tpval:{pval:.2e}",
            showarrow=False,


        )

        stayed_df = df[df["Accession"].isin(res_ids)]

        fig.add_trace(go.Scatter(
            x=df["log2FoldChange_RNA"],
            y=df["log2FoldChange_Protein"],
            mode="markers",
            text=df["hover"],
            hovertemplate=hovertemplate,
            customdata=cd
        ))
        cd = np.stack((stayed_df["padj_RNA"], stayed_df["padj_Protein"]),axis=1)

        fig.add_trace(go.Scatter(
            x=stayed_df["log2FoldChange_RNA"],
            y=stayed_df["log2FoldChange_Protein"],
            mode="markers",
            text=stayed_df["hover"],
            hovertemplate=hovertemplate,
            customdata=cd,
            name="Stayed at membrane"
        ))
        fig.update_yaxes(title="Protein l2FC")
        fig.update_xaxes(title="RNA l2FC")
        fig.write_html(output.html)


rule checkProteinLocalization:
    input:
        uniprot = rules.download_subcell.output.tsv,
        deseq_res = rules.joinTables.output.tsv,
    output:
        svg = "Pipeline/SubcellularLoc/{condition}vs{baseline}SubcellLoc.svg",
        tsv = "Pipeline/SubcellularLoc/{condition}vs{baseline}SubcellLoc.tsv",
        r_is_not_a_valid_programming_language_and_no_one_should_use_this_bullshit = "Pipeline/SubcellularLoc/{condition}vs{baseline}SubcellLocSameBullshitButCauseRCantReadDataframesItNeedsToBePreprocessed.tsv"

    script: "scripts/checkMembraneDistribution.py"



rule cpF3SubCellLoc:
    input:
        file = expand(rules.checkProteinLocalization.output.svg, condition=["M"], baseline=["C"])[0]
    output:
        svg = config["F3SubCellLoc"]["file"]
    shell:
        "cp {input.file} {output.svg}"


rule TransmembraneGSEA:
    input:
        tsv = rules.checkProteinLocalization.output.r_is_not_a_valid_programming_language_and_no_one_should_use_this_bullshit
    output:
        enriched="Pipeline/customGSEA/tsv/TransmembraneGSEA_p{puro}_c{condition}_vs_b{baseline}.tsv",
    params:
        terms=["Transmembrane", "NonTransmembrane"]
    script: "scripts/TransmembraneGSEA.R"


rule compareTotalCell:
    input:
        file = rules.tc_extractDESeqResult.output.result_table,
    output:
        svg="Pipeline/DFGAntrag/summary_total_cell_c_{condition}vs{baseline}.svg",
        html="Pipeline/DFGAntrag/summary_total_cell_c_{condition}vs{baseline}.html",
    run:
        from helper import compare_results
        cfg = config
        for (idx, (inputs, outputs)) in enumerate([(input.file, (output.svg, output.html))]):
            puro = "Total cell"
            fig = compare_results(
                [(puro, inputs)],
                sep="\t",
                condition=NAMEMapping[wildcards.condition],
                base=NAMEMapping[wildcards.baseline],
                padj_cutoff=config["primary"]["log2FCCutOff"],
                lfc_cutoff=config["primary"]["log2FCCutOff"],
                colors=(cfg["cytosolColor"], "#eaeaea", cfg["membraneColor"]),
                y="Name")
            fig.update_yaxes(title=None,showticklabels=False,showline=True,mirror=True,linewidth=1,linecolor="black")
            fig.update_xaxes(title="Fraction of transcriptome [%]",showline=True,mirror=True,linewidth=1,linecolor="black",gridcolor="black",dtick=20)
            print(fig.layout.legend)
            fig.update_layout(
                title=dict(text=f"{puro}: plus Puromycin vs minus Puromycin"),
                width=360,
                height=150,
                font=dict(size=10,color="black"),
                margin=dict(r=10,l=10,t=50,b=20),
                bargap=0.3,
                plot_bgcolor='white',
                legend=dict(
                    x=1,# Adjust the x position of the legend
                    y=0.5,# Adjust the y position of the legend
                    xanchor='right',# Anchor the legend to the right side
                    yanchor='middle',# Anchor the legend to the top
                    xref="container",
                    yref="container"
                )
            )
            fig.update_layout(
                yaxis=dict(domain=[0.6, 1]),# Adjust the y-axis domain
            )
            fig.write_image(outputs[0])
            fig.write_html(outputs[1])



rule compareResults:
    input:
        file = rules.other_extractDESeqResult.output.result_table,
        file2 = rules.second_extractDESeqResult.output.result_table,

    output:
        svg = "Pipeline/DFGAntrag/summary_minus_c_{condition}vs{baseline}.svg",
        html = "Pipeline/DFGAntrag/summary_minus_c_{condition}vs{baseline}.html",
        svgplus = "Pipeline/DFGAntrag/summary_plus_c_{condition}vs{baseline}.svg",
        htmlplus = "Pipeline/DFGAntrag/summary_plus_c_{condition}vs{baseline}.html",
    run:
        from helper import compare_results
        for (idx, (inputs, outputs)) in enumerate(((input.file, (output.svg, output.html)), (input.file2, (output.svgplus, output.htmlplus)))):
            puro = "Puromycin-" if idx == 0 else "Puromycin+"
            fig = compare_results(
                [(puro, inputs)],
                sep="\t",
                condition=NAMEMapping[wildcards.condition],
                base=NAMEMapping[wildcards.baseline],
                padj_cutoff=config["primary"]["log2FCCutOff"],
                lfc_cutoff=config["primary"]["log2FCCutOff"],
                colors=("#002695", "#eaeaea", "#2ca02c"),
                y = "Name")
            fig.update_yaxes(title=None, showticklabels=False, showline=True, mirror=True, linewidth=1, linecolor="black")
            fig.update_xaxes(title="Fraction of transcriptome [%]", showline=True, mirror=True, linewidth=1, linecolor="black", gridcolor="black", dtick=20)
            fig.update_layout(
                title=dict(text=f"{puro}: C vs M"),
                width=624/2,
                height=150,
                font=dict(size=10, color="black"),
                margin=dict(r=10, l=10, t=50, b=20),
                bargap=0.3,
                plot_bgcolor='white',
                legend=dict(
                    x=1,# Adjust the x position of the legend
                    y=0.5,# Adjust the y position of the legend
                    xanchor='right',# Anchor the legend to the right side
                    yanchor='middle',# Anchor the legend to the top
                    xref="container",
                    yref="container"
                )
            )
            fig.update_layout(
                yaxis=dict(domain=[0.6, 1]),# Adjust the y-axis domain
            )
            fig.write_image(outputs[0])
            fig.write_html(outputs[1])

rule runOnLeastDiffRegGenes:
    input:
        files = expand(rules.other_extractDESeqResult.output.result_table,  zip, condition=config["primary"]["conditions"], baseline=config["primary"]["baselines"])
    output:
        file =   os.path.join(config["secondary"]["RUN_DIR"], "PipelineData/Multilevel/intersectingGenes.tsv")
    run:
        import pandas as pd
        sets = []
        for file in input.files:
            df = pd.read_csv(file, sep="\t")
            df = df[(df["log2FoldChange"] <= config["multilevel"]["lfc_cutoff"]) & (df["log2FoldChange"] >= -config["multilevel"]["lfc_cutoff"]) ]
            df = df[df["padj"] >= config["multilevel"]["padj_cutoff"]]
            sets.append(set(df.index))
        hs_genes = config["multilevel"]["additional_housekeeping_genes"]
        if len(hs_genes) > 0:
            for gene in hs_genes:
                assert gene in df.index, f"{gene} not found in index"
            sets.append(hs_genes)
        if config["multilevel"]["mode"] == "intersect":
            sets = set.intersection(*sets)
        else:
            raise NotImplementedError("other mode not implemented yet")
        df = df[df.index.isin(sets)]
        df.to_csv(output.file, sep="\t")

#config["lightDarkLight"]["housekeeping"] = str(rules.runOnLeastDiffRegGenes.output.file)
#config["lightDarkLight"]["use-housekeeping"] = True


#config["lightDarkDark"]["housekeeping"] = str(rules.runOnLeastDiffRegGenes.output.file)
#config["lightDarkDark"]["use-housekeeping"] = True

module ld_light_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkLight"]

use rule * from ld_light_workflow as lightDarkLight_*

module ld_dark_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkDark"]

use rule * from ld_dark_workflow as lightDarkDark_*



module ld_light_puro_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkLightPuro"]

use rule * from ld_light_puro_workflow as lightDarkLightPuro_*

module ld_dark_puro_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkDarkPuro"]

use rule * from ld_dark_puro_workflow as lightDarkDarkPuro_*

module ld_light_tc_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkLightTC"]

use rule * from ld_light_tc_workflow as lightDarkLightTC_*

module ld_light_vs_light_tc_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkLightvsLightTC"]

use rule * from ld_light_vs_light_tc_workflow as lightDarkLightvsLightTC_*

module ld_light_vs_light_plus_puro_tc_workflow:
    snakefile: snakefile
    config: config["RUNS"]["lightDarkLightvsLightPlusPuroTC"]

use rule * from ld_light_vs_light_plus_puro_tc_workflow as lightDarkLightvsLightPlusPuroTC_*

rule extractOverlap:
    input:
        deseq1 = rules.other_extractDESeqResult.output.result_table,
        deseq2 = rules.lightDarkLight_extractDESeqResult.output.result_table,
        deseq3 = rules.lightDarkDark_extractDESeqResult.output.result_table,
        lp = rules.second_extractDESeqResult.output.result_table,
        ldlp = rules.lightDarkLightPuro_extractDESeqResult.output.result_table,
        lddp = rules.lightDarkDarkPuro_extractDESeqResult.output.result_table,
    output:
        svg = "Pipeline/overlap/c{condition}_vs_b{baseline}_overlap.svg",
        html = "Pipeline/overlap/c{condition}_vs_b{baseline}_overlap.html"
    conda: "envs/matplotlibPlotly.yml"
    script:
        "scripts/getVennDiagram.py"


rule customGSEA:
    input:
        deseq_results = "Pipeline/RUNS/{run_dir}/PipelineData/DESeqResults/DESeqResult_c{condition}_vs_b{baseline}.tsv"
    output:
        enriched="Pipeline/customGSEA/tsv/rd{run_dir}_c{condition}_vs_b{baseline}_UTREnrich.tsv",
        gsdata="Pipeline/customGSEA/tsv/rd{run_dir}_c{condition}_vs_b{baseline}_GSDATA.tsv",
    params:
        Rlib = config["primary"]["Rlib"]
    conda: "envs/clusterprofiler.yml"
    script:
        "scripts/GSEA.R"

rule plotGSEA:
    input:
        df = rules.customGSEA.output.gsdata,
        info = rules.customGSEA.output.enriched
    output:
        html = "Pipeline/customGSEA/plot/rd{run_dir}_c{condition}_vs_b{baseline}_UTREnrich.html",
    conda: "envs/matplotlibPlotly.yml"
    script:
        "scripts/GSEAPlot.py"



rule extractFasta:
    input:
        gff = config["gff"],
        fasta = config["fasta"]
    output:
        fasta = "Pipeline/sequences.fa",
    conda: "envs/bedtools.yml"
    shell:
        """
        bedtools getfasta -s -fi {input.fasta} -bed {input.gff} > {output.fasta}
        """


rule plotGSEARibosome:
    input:
        tsv = expand(rules.other_GSEAGO.output.enriched, condition=["M"], baseline=["C"])[0],
        info_tsv = expand(rules.other_GSEAGO.output.gsdata, condition=["M"], baseline=["C"])[0]
    params:
        base_name="Cytosol",
        condition_name="Membrane",
        layout=config["F3GSEA"]["layout"]
    output:
        svg = config["F3GSEA"]["file"]
    conda: "envs/matplotlibPlotly.yml"
    script:
        "scripts/GOGSEA.py"

rule addFastatoTSV:
    input:
        fasta = rules.extractFasta.output.fasta,
        gff = config["gff"],
        tsv = expand(rules.checkProteinLocalization.output.tsv, condition="M", baseline="C")[0]
    output:
        tsv="Pipeline/TableWithSequences.tsv",
    conda: "envs/biopython.yml"
    script:
        "scripts/addFastaToTSV.py"



rule getFastaForMotif:
    input:
        tsv = rules.addFastatoTSV.output.tsv
    output:
        bg = "Pipeline/stayedMotif/bg_seqs_t{type}.fa",
        stayed = "Pipeline/stayedMotif/stayed_seqs_t{type}.fa"
    conda: "envs/biopython.yml"
    script: "scripts/extractFastaFromTable.py"


rule memeMotifStayed:
    input:
        stayed = rules.getFastaForMotif.output.stayed,
        bg = rules.getFastaForMotif.output.bg
    output:
        dir = directory("Pipeline/stayedMotif/MEME_t{type}/"),
    conda:
        "envs/meme.yml"
    shell:
        """
        streme --p {input.stayed} --o {output.dir}  --evalue --minw 3 --align center --rna --nmotifs 5 --maxw 6
        """



rule checkMotifDistribution:
    input:
        tsv = rules.addFastatoTSV.output.tsv,
    conda: "envs/biopython.yml"

    output:
        tsv = "Pipeline/Uracil/UracilDistribution.tsv"
    script: "scripts/checkMotifDistribution.py"

#
# rule prepareForGSEAMotif:
#     input:
#         motifs = rules.checkMotifDistribution.output.tsv,
#         deseq = expand(rules.joinTables.output.output_for_R_because_this_stupid_language_cannot_even_read_tsv_files, condition=["M"], baseline=["C"])
#     output:
#         file = "Pipeline/MotifGSEA/preparedFile_t{type}.tsv"
#     run:
#         import pandas as pd
#         import os
#         deseq = pd.read_csv(input.deseq[0], sep="\t")
#         motifs = pd.read_csv(input.motifs, sep="\t", comment="#")
#
#         for motif_id in motifs.motif_ID.unique():
#             seqs = motifs[motifs["motif_ID"] == motif_id]["seq_ID"]
#             deseq["Motif" + motif_id.replace("-", "_") ] = deseq.locus_tag.isin(seqs)
#         deseq = deseq[deseq["mRNA"]]
#         deseq.to_csv(output.file, sep="\t", index=False)
#
#
#
# rule motifGSEA:
#     input:
#         tsv = rules.prepareForGSEAMotif.output.file
#     output:
#         enrichedFile = "Pipeline/MotifGSEA/enrichedFile_t{type}.tsv",
#         gsFile = "Pipeline/MotifGSEA/gsFile_t{type}.tsv"
#     conda: "envs/clusterprofiler.yml",
#     script: "scripts/motifGSEA.R"

# rule plotMotifGSEA:
#     input:
#         df = rules.motifGSEA.output.gsFile,
#         info = rules.motifGSEA.output.enrichedFile
#     output:
#         html = "Pipeline/MotifGSEA/GSEAMotif_t{type}.html",
#     conda: "envs/matplotlibPlotly.yml"
#     script:
#         "scripts/GSEAPlot.py"


rule ntFreqs:
    input:
        stayed=rules.getFastaForMotif.output.stayed,
        bg=rules.getFastaForMotif.output.bg
    output:
        dir = "Pipeline/stayedMotif/NtFreqs_t{type}.tsv",
    conda: "envs/biopython.yml"
    script:
        "scripts/calcNTFreq.py"





rule addRBP3Data:
    input:
        rbp3  = "initialData/rbp3data.tsv",
        joined = rules.addFastatoTSV.output.tsv
    output:
        trash_tsv_for_R = "Pipeline/RBP3Check/trashTSVForR.tsv",
        useful_tsv = "Pipeline/RBP3Check/rbp3joined.tsv",
        html = "Pipeline/RBP3Check/rbp3joined.html",
    script: "scripts/checkRBP3Data.py"

rule rnatypeGSEA:
    input:
        tsv = expand(rules.joinTables.output.output_for_R_because_this_stupid_language_cannot_even_read_tsv_files, condition="M", baseline="C")[0]
    output:
        enrichedFile = "Pipeline/rnaTypeGSEA/gseaEnrich_puro{puro}.tsv",
        gsFile = "Pipeline/rnaTypeGSEA/gseaInfo_puro{puro}.tsv"
    params:
        terms = config["RNA_TYPE_COLORS"].keys()
    conda: "envs/clusterprofiler.yml",
    script: "scripts/rnatypeGSEA.R"



rule rbp3GSEA:
    input:
        tsv = rules.addRBP3Data.output.trash_tsv_for_R,
    output:
        enrichedFile = "Pipeline/RBP3Check/gseaEnrich_puro{puro}.tsv",
        gsFile = "Pipeline/RBP3Check/gseaInfo_puro{puro}.tsv"
    conda: "envs/clusterprofiler.yml",
    script: "scripts/RBP3GSEA.R"

rule plotRBP3GSEA:
    input:
        df=rules.rbp3GSEA.output.gsFile,
        info=rules.rbp3GSEA.output.enrichedFile,
    output:
        html = "Pipeline/RBP3Check/RBP3_GSEA_{puro}.html",
        svg = "Pipeline/RBP3Check/RBP3_GSEA_{puro}.svg",
    params:
        figureConfig="RBP3GSEA",
        show_zero_lfc=True,
        condition_name=None,
        base_name=None
    conda: "envs/matplotlibPlotly.yml"
    script:
        "scripts/RBP3GSEAPlot.py"

rule plotRBP3Venns:
    input:
        tsv = rules.addRBP3Data.output.useful_tsv
    output:
        html = "Pipeline/RBP3Check/RBP3_Venn.html",
        svg = "Pipeline/RBP3Check/RBP3_Venn.svg",
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/rbp3Venns.py"


rule joinRBP3Plot:
    input:
        f1 = rules.plotRBP3Venns.output.svg,
        f2 = expand(rules.plotRBP3GSEA.output.svg, puro=["plusPuromycin"])[0],
    output:
        svg = "Pipeline/Paper/RBP3Figure.svg"
    conda: "envs/svgutils.yml"
    script: "scripts/joinPlots.py"


rule extractRBP3Seqs:
    input:
        tsv = rules.addRBP3Data.output.useful_tsv
    params:
        type= "RBP3_Target"
    output:
        fasta = "Pipeline/RBP3Check/rbp3Seqs.fa"
    conda: "envs/biopython.yml"
    script: "scripts/extractRBP3Fasta.py"

rule fimoRBP3:
    input:
        fasta = rules.extractRBP3Seqs.output.fasta,
        motifs = rules.memeMotifStayed.output.dir
    output:
        dir = directory("Pipeline/RBP3Check/FIMO_t{type}/")
    conda:
        "envs/meme.yml"
    shell:
        "fimo  --oc {output.dir}  --thresh 1.0E-4 --bfile --motif-- {input.motifs}/streme.txt {input.fasta}"



rule minusPVolcanoForFig1:
    input:
        tsv = expand(rules.joinTables.output.tsv, condition="M", baseline="C")[0],
        gsea = expand(rules.rnatypeGSEA.output.enrichedFile, puro="minusPuromycin")[0]
    params:
        condition="Membrane",
        baseline= "Cytosol",
        figureConfig="F1Volcano"
    output:
        svg = config["F1Volcano"]["file"],
        html = "Pipeline/Paper/Subfigures/F1Volcano.html"
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/refineVolcano.py"

rule refinePercentagePlot:
    input:
        tsv = expand(rules.other_extractDESeqResult.output.result_table,condition=["M"], baseline=["C"])[0]
    output:
        svg = config["F1Percent"]["file"]
    params:
        condition="Membrane",
        baseline="Cytosol",
        figureConfig="F1Percent",
        name="FOO",
        log2FCCutOff=config["primary"]["log2FCCutOff"],
        pAdjCutOff=config["primary"]["pAdjCutOff"]

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/percentagePlot.py"

rule PercentagePlotF4:
    input:
        tsv = expand(rules.tc_extractDESeqResult.output.result_table,condition=["minusP"], baseline=["plusP"])[0]
    output:
        svg = config["F4Percent"]["file"],
    params:
        condition="- Puromycin",
        baseline="+ Puromycin",
        figureConfig="F4Percent",
        name="FOO",
        log2FCCutOff=config["primary"]["log2FCCutOff"],
        pAdjCutOff=config["primary"]["pAdjCutOff"]

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/percentagePlot.py"


rule PercentagePlotF5:
    input:
        tsv = expand(rules.second_extractDESeqResult.output.result_table,condition=["M"], baseline=["C"])[0]
    output:
        svg = config["F5Percent"]["file"]
    params:
        condition="Membrane",
        baseline="Cytosol",
        figureConfig="F1Percent",
        name="FOO",
        log2FCCutOff=config["primary"]["log2FCCutOff"],
        pAdjCutOff=config["primary"]["pAdjCutOff"]

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/percentagePlot.py"

rule F1GSEA:
    input:
        df=rules.rnatypeGSEA.output.gsFile,
        info=rules.rnatypeGSEA.output.enrichedFile,
    output:
        html = "Pipeline/rnaTypeGSEA/rnaTypeGSEA_{puro}.html",
        svg="Pipeline/Paper/Subfigures/rnaTypeGSEA_{puro}.svg",
    conda: "envs/matplotlibPlotly.yml"
    params:
        figureConfig="F1GSEA",
        show_zero_lfc=True,
        condition_name=None,
        base_name=None
    script:
        "scripts/RBP3GSEAPlot.py"

rule txtToPathWorkflow:
    input:
        file = "Workflow.svg"
    output:
        svg = config["F1Workflow"]["file"]
    shell: "inkscape --export-text-to-path -o {output.svg} {input.file}"


rule joinFig1:
    input:
        volcano = rules.minusPVolcanoForFig1.output.svg,
        summary = rules.refinePercentagePlot.output.svg,
        workflow = rules.txtToPathWorkflow.output.svg,
        gsea = expand(rules.F1GSEA.output.svg, puro="minusPuromycin")[0]
    params:
        figures = config["F1"],
        ges_y = 520
    conda: "envs/svgutils.yml"
    output:
        svg = "Pipeline/Paper/Figure1/Figure1.svg"
    script:
        "scripts/paperFigure.py"

rule Figure1toPNG:
    input:
        svg = rules.joinFig1.output.svg
    output:
        png = "Pipeline/Paper/Figure1/Figure1.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"






rule plotCoveragesFig2:
    input:
        design = "initialData/testDesign.tsv",
        gff = "initialData/annotationFiltered.gff",
        sf = rules.other_runDESeq.output.size_factors,
        db= "Pipeline/GFFDB/gff_database.db",
    output:
        svg = config["MinusPuroCoverage"]["file"]
    conda: "envs/matplotlibPlotly.yml"

    script: "scripts/plotCoveragesFigure2.py"




# rule plotGOTermEnrichments:
#     input:
#         stayed = expand(rules.GOTermEnrichmentStayed.output, condition=["M"], baseline=["C"]),
#     output:
#         svg =


rule processFISHImages:
    output:
        tsv = "Pipeline/FISH/all_distributions.tsv",
        spots_df = "Pipeline/FISH/single_spots.tsv"
    threads: 10
    notebook:
        "scripts/generateFISHDistribution.ipynb"


rule plotFISHFigure2:
    input:
        tsv = rules.processFISHImages.output.tsv
    output:
        svg = config["F2FISH"]["file"]
    params:
        probes = [("psaA", "untreated", "Replicate 1"), ("cupA", "untreated", "Replicate 1"), ("rnpB", "untreated", "Replicate 1"), ("fusB", "Puromycin", "Replicate 1") ],
        layout= config["F2FISH"]["layout"],
        subplot=dict(vertical_spacing=0.05, column_titles=["psaA", "cupA", "rnpB", "fusB"])

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotFISH.py"


rule joinFig2:
    input:
        coverage = rules.plotCoveragesFig2.output.svg,
        fish = rules.plotFISHFigure2.output.svg
    params:
        figures = config["F2"],
        ges_y = 450
    conda: "envs/svgutils.yml"
    output:
        svg = "Pipeline/Paper/Figure2/Figure2.svg"
    script:
        "scripts/paperFigure.py"

rule Figure2toPNG:
    input:
        svg = rules.joinFig2.output.svg
    output:
        png = "Pipeline/Paper/Figure2/Figure2.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"


rule plotFISHFigure3:
    input:
        tsv=rules.processFISHImages.output.tsv
    output:
        svg = config["F3FISH"]["file"]
    params:
        probes =[("16S rRNA", "untreated", "Replicate 1")],
        layout= config["F3FISH"]["layout"],
        subplot=dict(vertical_spacing=0.05, horizontal_spacing=0.05,)

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotFISHFigure3.py"


# rule plotFISHFigure4OLD:
#     input:
#         tsv=rules.processFISHImages.output.tsv
#     output:
#         svg = config["F4FISH"]["file"]
#     params:
#         probes =[("16S rRNA", "untreated", "Replicate 1")],
#         layout= config["F4FISH"]["layout"],
#         subplot=dict(vertical_spacing=0.05, horizontal_spacing=0.05,)
#
#     conda: "envs/matplotlibPlotly.yml"
#     script: "scripts/OLDplotFISHFigure4.py"

rule F4ClusteredGO:
    input:
        tsv = expand(rules.tc_ClusterSemSim.output.table, condition="plusP", baseline="minusP", updown="up")[0],
    params:
        layout = config["F4ClusteredGO"]["layout"]
    output:
        svg = config["F4ClusteredGO"]["file"]
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotClusteredEnrichment.py"

rule plotQPCR:
    input:
        tsv = "initialData/ActualQPCRData.csv"
    conda: "envs/matplotlibPlotly.yml"
    output:
        svg = config["F4QTPCR"]["file"]
    script:
        "scripts/QTPCR.py"

rule plotF4cupAIntensities:
    input:
        tsv=rules.processFISHImages.output.spots_df,
    output:
        svg=config["F4CupAIntensities"]["file"],
        tsv= "Pipeline/Paper/Figure4/IntensitiesTestPval.tsv"

    script:
        "scripts/plotFISHIntensity.py"




rule plotVennsFigure5:
    input:
        tsv=expand(rules.joinTables.output.tsv,condition="M",baseline="C")[0],
    conda: "envs/matplotlibPlotly.yml"
    output:
        svg = config["Fig5Venn"]["file"],
    script: "scripts/VennOverlap.py"



rule plotVolcanoFigure5:
    input:
        tsv=expand(rules.second_extractDESeqResult.output.result_table,condition="M",baseline="C")[0],
        go=expand(rules.second_GOEnrichment.output.up, condition="M", baseline="C")[0],
        add_data = "initialData/locusTagtoGeneName.csv"
    output:
        svg = config["VolcanoFigure5"]["file"],
        html="Pipeline/Paper/Figure5/F5Volcano.html"
    params:
        figureConfig = "VolcanoFigure5"
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/VolcanoFigure5.py"






rule plotFISHFigure4:
    input:
        tsv=rules.processFISHImages.output.tsv
    output:
        svg = config["F4FISH"]["file"]
    params:
        probes =[("psbA2","untreated","Replicate 1"), ("psbA2","Puromycin","Replicate 1"), ("16S rRNA", "untreated", "Replicate 1"), ("16S rRNA", "Puromycin", "Replicate 1")],
        layout= config["F4FISH"]["layout"],
        subplot=dict(vertical_spacing=0.05, horizontal_spacing=0.05,
                              column_titles=["- Puromycin", "+ Puromycin", "- Puromycin", "+ Puromycin"])

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotFISHFigure4.py"


rule plotFISHFigure5:
    input:
        tsv=rules.processFISHImages.output.tsv,
        intensities=rules.processFISHImages.output.spots_df
    output:
        svg = config["F5FISH"]["file"]
    params:
        probes =[("fusB","untreated","Replicate 1"), ("fusB","Puromycin","Replicate 1"), ("psaA", "untreated", "Replicate 1"), ("psaA", "Puromycin", "Replicate 1")],
        layout= config["F5FISH"]["layout"],
        subplot=dict(vertical_spacing=0.05, horizontal_spacing=0.05,
                              column_titles=["- Puromycin", "+ Puromycin", "- Puromycin", "+ Puromycin"])

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotFISHFigure4.py"




rule plotFISHFigure6:
    input:
        tsv=rules.processFISHImages.output.tsv,
    output:
        svg = config["F6FISH"]["file"]
    params:
        probes =[("psbA2","untreated","Replicate 1"), ("psbA2","Puromycin","Replicate 1"), ("rpl10_12", "untreated", "Replicate 1"), ("rpl10_12", "Puromycin", "Replicate 1")],
        layout= config["F6FISH"]["layout"],
        subplot=dict(vertical_spacing=0.05, horizontal_spacing=0.05,
                              column_titles=["- Puromycin", "+ Puromycin", "- Puromycin", "+ Puromycin"])

    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotFISHFigure6.py"

rule plotCoveragesFig6:
    input:
        design = "initialData/testDesign.tsv",
        gff = "initialData/annotationFiltered.gff",
        sf = rules.third_runDESeq.output.size_factors,
    output:
        svg = config["F6Coverage"]["file"],
        db = "Pipeline/GFFDB/gff_database.db",
    conda: "envs/matplotlibPlotly.yml"

    script: "scripts/plotCoveragesFigure6.py"

rule plotGOFigure6:
    input:
        tsv=expand(rules.second_ClusterSemSim.output.table,condition="M",baseline="C",updown="up")[0],
    params:
        layout=config["F6ClusteredGO"]["layout"]
    output:
        svg=config["F6ClusteredGO"]["file"]
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotClusteredEnrichment.py"







rule joinFig6:
    input:
        files = [config[val]["file"] for key, val in config["F6"].items()]
    params:
        figures = config["F6"],
        ges_y = 750
    conda: "envs/svgutils.yml"
    output:
        svg = "Pipeline/Paper/Figure6/Figure6.svg"
    script:
        "scripts/paperFigure.py"


rule Figure6toPNG:
    input:
        svg = rules.joinFig6.output.svg
    output:
        png = "Pipeline/Paper/Figure6/Figure6.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"

rule joinFig5:
    input:
        files=[config[val]["file"] for key, val in config["F5"].items()]
    params:
        figures = config["F5"],
        ges_y = 700
    conda: "envs/svgutils.yml"
    output:
        svg = "Pipeline/Paper/Figure5/Figure5.svg"
    script:
        "scripts/paperFigure.py"

rule checkTranscriptAbundanceGO:
    input:
        go  =expand(rules.other_GOEnrichment.output.up, condition="M", baseline="C")[0],
        tsv =expand(rules.tc_extractDESeqResult.output.result_table,condition="plusP",baseline="minusP")[0],
    output:
        tsv = "Pipeline/SHIT/abundanceOfPlasmaMembrane.tsv"
    script: "scripts/checkTranscriptAbundanceGO.py"

rule Figure5toPNG:
    input:
        svg = rules.joinFig5.output.svg
    output:
        png = "Pipeline/Paper/Figure5/Figure5.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"

rule plotVolcanoFigure4:
    input:
        tsv=expand(rules.tc_extractDESeqResult.output.result_table,condition="plusP",baseline="minusP")[0],
        go=expand(rules.tc_GOEnrichment.output.up, condition="plusP", baseline="minusP")[0],
        add_data= "initialData/locusTagtoGeneName.csv"
    params:
        figureConfig = "F4Volcano"
    output:
        svg = config["F4Volcano"]["file"],
        html="Pipeline/Paper/Figure4/F4Volcano.html"
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/VolcanoFigure5.py"

rule plotGOTermEnrichmentFig4:
    input:
        tsv =expand(rules.tc_GOEnrichment.output.up, condition="plusP", baseline="minusP")[0],
    output:
        svg = "Pipeline/Paper/Figure4/F4Enrichment.svg",
        html = "Pipeline/Paper/Figure4/F4Enrichment.html",
    conda: "envs/matplotlibPlotly.yml"

    script: "scripts/enrichmentTotalCell.py"

rule plotFULLGOTermEnrichment:
    input:
        tsv1 =expand(rules.other_GOEnrichment.output.up, condition="M", baseline="C")[0],
        tsv2= expand(rules.second_GOEnrichment.output,condition=["M"],baseline=["C"])[0],
        tsv3= expand(rules.GOTermEnrichmentStayed.output,condition=["M"],baseline=["C"])[0],
    output:
        svg = "Pipeline/Paper/Supplementary/FullEnrichment.svg",
        html = "Pipeline/Paper/Supplementary/FullEnrichment.html",
    conda: "envs/matplotlibPlotly.yml"

    script: "scripts/plotEnrichmentPlots.py"


rule F3ClusteredGO:
    input:
        tsv = expand(rules.other_ClusterSemSim.output.table, condition="M", baseline="C", updown="up")[0],
    params:
        layout = config["F3ClusteredGO"]["layout"]
    output:
        svg = config["F3ClusteredGO"]["file"]
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotClusteredEnrichment.py"

rule F3QTPCR:
    output:
        svg = "Pipeline/Paper/Figure3/F3QTPCR.svg",
    script: "scripts/QTPCR.py"


rule joinFig3:
    input:
        files = [config[val]["file"] for key, val in config["F3"].items()]
    params:
        figures = config["F3"],
        ges_y = 500
    conda: "envs/svgutils.yml"
    output:
        svg = "Pipeline/Paper/Figure3/Figure3.svg"
    script:
        "scripts/paperFigure.py"

rule Figure3toPNG:
    input:
        svg = rules.joinFig3.output.svg
    output:
        png =  "Pipeline/Paper/Figure3/Figure3.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"


rule joinFig4:
    input:
        files=[config[val]["file"] for key, val in config["F4"].items()]
    params:
        figures = config["F4"],
        ges_y = 800
    conda: "envs/svgutils.yml"
    output:
        svg = "Pipeline/Paper/Figure4/Figure4.svg"
    script:
        "scripts/paperFigure.py"

rule Figure4toPNG:
    input:
        svg = rules.joinFig4.output.svg
    output:
        png =  "Pipeline/Paper/Figure4/Figure4.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"



rule all:
    default_target: True
    input:
        all = rules.other_all.input,
        all2 = rules.second_all.input,
        tc_alls = rules.tc_all.input,
        motif = rules.meme_allMotif.input,
        motif2 = rules.meme2_allMotif.input,
        motifldl = rules.meme_ldl_allMotif.input,
        motifldd = rules.meme_ldd_allMotif.input,
        motifldlp = rules.meme_ldlpuro_allMotif.input,
        motiflddp = rules.meme_lddpuro_allMotif.input,
        tu = expand(rules.createTUTable.output, condition=["M"], baseline=["C"]),
        tu2 = expand(rules.createTUTable2.output, condition=["M"], baseline=["C"]),
        #tu3 = expand(rules.plotHistogramm.output,zip, condition=["M"], baseline=["C"]),
        stayedGO = expand(rules.GOTermEnrichmentStayed.output, condition=["M"], baseline=["C"]),
        transintramembrane = expand(rules.joinTables.output, condition=["M"], baseline=["C"]),
        dfr = expand(rules.checkProteinLocalization.output, condition=["M"], baseline=["C"]),
        size_factors = rules.third_runDESeq.output,
        #protloc = expand(rules.correlationRNAvsProteins.output, condition="M", baseline="C"),
        summary_fig = expand(rules.compareResults.output, condition="M", baseline="C"),
        summary_fig2 = expand(rules.compareTotalCell.output, condition="plusP", baseline="minusP"),
        #overlap= expand(rules.extractOverlap.output, zip, condition=config["primary"]["conditions"],baseline=config["primary"]["baselines"]),
        ldl = rules.lightDarkLight_all.input,
        ldd = rules.lightDarkDark_all.input,
        ldl_puro = rules.lightDarkLightPuro_all.input,
        ldd_puro = rules.lightDarkDarkPuro_all.input,
        ldd_tc = rules.lightDarkLightTC_all.input,
        ldd_vsl_tc = rules.lightDarkLightvsLightTC_all.input,
        ldd_vsl_plusp_tc = rules.lightDarkLightvsLightPlusPuroTC_all.input,
        test = expand(rules.plotGSEA.output, condition="M", baseline="C", run_dir=[os.path.basename(config["RUNS"][run]["RUN_DIR"]) for run in config["RUNS"] if "TC" not in run]),
        subcell = rules.download_subcell.output,
        rbp3 = expand(rules.plotRBP3GSEA.output, puro=["plusPuromycin", ]),
        rbp3venn = rules.joinRBP3Plot.output,
        rbp3fimo = expand(rules.fimoRBP3.output, type=["Cellular_thylakoid_membrane", "Cell_membrane"]),
        #prdwep = expand(rules.checkMotifDistribution.output, type=["Cellular_thylakoid_membrane", "Cell_membrane"]),
        fig1 = rules.joinFig1.output,
        fig1png = rules.Figure1toPNG.output,
        rnatypegsea = expand(rules.rnatypeGSEA.output, puro=["minusPuromycin", ]),
        cov2 = rules.plotCoveragesFig2.output,
        f5volc = rules.plotVolcanoFigure5.output,
        f4volc = rules.plotVolcanoFigure4.output,
        f4enrich = rules.plotGOTermEnrichmentFig4.output,
        f6enrich = rules.plotFULLGOTermEnrichment.output,
        f5venn = rules.plotVennsFigure5.output,
        fish = rules.Figure2toPNG.output,
        fig5 = rules.Figure5toPNG.output,
        fish3 = rules.Figure6toPNG.output,
        f3 = rules.Figure3toPNG.output,
        t1 = rules.Table1.output,
        f3ribo = rules.plotGSEARibosome.output,
        qpcr = rules.Figure4toPNG.output,
        clusteredF4 = rules.F4ClusteredGO.output,
        summary = expand(rules.summarizeResult.output, condition="M", baseline="C"),
        transmemgsea = expand(rules.TransmembraneGSEA.output.enriched, puro=["minusPuromycin", "plusPuromycin"], condition="M", baseline="C"),
        f5perc = config["F5Percent"]["file"],
        f5fish = config["F5FISH"]["file"],
        #tpmemabundace = rules.checkTranscriptAbundanceGO.output


