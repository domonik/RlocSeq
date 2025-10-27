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

with open(config["why-drop-replicate-file"], "r") as handle:
    config["why_drop_replicate"] = yaml.safe_load(handle)



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

module why_drop_workflow:
    snakefile: snakefile
    config: config["why_drop_replicate"]

use rule * from why_drop_workflow as why_drop_*


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


rule joinTables:
    input:
        deseqResMinus = expand(rules.other_extractDESeqResult.output.result_table, condition="M", baseline="C")[0],
        deseqResPlus = expand(rules.second_extractDESeqResult.output.result_table, condition="M", baseline="C",)[0],
        deseqResTC = expand(rules.tc_extractDESeqResult.output.result_table, condition="plusP", baseline="minusP")[0],
        add_data = "initialData/locusTagtoGeneName.csv",
        tu_info = rules.coLocateTranscripts.output.tus,
        uniprot= rules.download_subcell.output.tsv,

    output:
        tsv = "Pipeline/Summary/joinedDESeq_Table.tsv",
        output_for_R_because_this_stupid_language_cannot_even_read_tsv_files = "Pipeline/Summary/joinedDESeqForstupidR.tsv",
        stayed = "Pipeline/Summary/stayedAtMembrane.tsv",
    run:
        df_minus = pd.read_csv(input.deseqResMinus, sep="\t")
        df_plus = pd.read_csv(input.deseqResPlus, sep="\t")
        df_puro = pd.read_csv(input.deseqResTC, sep="\t")
        df = df_minus.join(df_plus, rsuffix="_plusPuromycin", lsuffix="_minusPuromycin", how="outer")
        df_puro.columns = [f"{col}_tc" for col in df_puro.columns]
        df = df.join(df_puro, rsuffix="_totalCell", lsuffix="_totalCell", how="outer")
        df = df[~df.index.str.contains("6803r")]

        # Here uniprot is joined:
        x = ("Cytoplasm", "Membrane", "Cell membrane", "Cell inner membrane", "Cellular thylakoid membrane")

        uniprot = pd.read_csv(input.uniprot,sep="\t")
        subcell = "Subcellular location [CC]"
        uni_col = "Gene Names (ordered locus)"
        uniprot = uniprot.rename({uni_col: "locus_tag"},axis=1)
        uniprot["locus_tag"] = uniprot["locus_tag"].str.split("; ")
        uniprot = uniprot.explode('locus_tag').reset_index(drop=True)

        uniprot = uniprot[["locus_tag", "Gene Names", subcell, "Transmembrane", "Sequence"]]  # dont need intramembrane since there are zero annotated
        uniprot['Cell membrane'] = uniprot[subcell].str.contains("Cell membrane")
        uniprot["Cell inner membrane"] = uniprot[subcell].str.contains("Cell inner membrane")
        uniprot['Cellular thylakoid membrane'] = uniprot[subcell].str.contains("Cellular thylakoid membrane")
        uniprot['Cytoplasm'] = (uniprot[subcell].str.contains("Cytoplasm") == True)
        uniprot['Cytoplasmic side'] = (uniprot[subcell].str.contains("Cytoplasmic side") == True)
        uniprot.loc[uniprot[subcell].isna(), "Cytoplasm"] = pd.NA
        uniprot['Membrane'] = (uniprot[subcell].str.contains("Membrane") == True) | uniprot['Cellular thylakoid membrane'] | uniprot[
            'Cell membrane'] | uniprot["Cell inner membrane"] | uniprot["Cytoplasmic side"]
        uniprot.loc[uniprot[subcell].isna(), "Membrane"] = pd.NA

        df["merge_tag"] = df.index
        df.loc[df["merge_tag"].str.contains("UTR"), "merge_tag"] = df.loc[df[
            "merge_tag"].str.contains("UTR"), "merge_tag"].str[:-5]
        old_index = df.index

        df = df.merge(uniprot, left_on="merge_tag", right_on="locus_tag", how="left")
        df.index = old_index

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
        tsv = rules.joinTables.output.tsv,
    output:
        csv = "Pipeline/Paper/Table1/Table1.csv"
    run:
        df = pd.read_csv(input.tsv, sep="\t")
        df["TransmembraneDomain"] = (~pd.isna(df["Transmembrane"])) == True
        names = ["mRNA", "5UTR", "3UTR", "tRNA", "asRNA", "Photosystem mRNA", "Ribosome mRNA", "CRISPR", "ncRNA", "TransmembraneDomain"]
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
        up = "Pipeline/Summary/GO/StayedAtMembrane_ud_up.tsv",
        down = "Pipeline/Summary/GO/StayedAtMembrane_ud_down.tsv"
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



rule checkProteinLocalization:
    input:
        deseq_res = rules.joinTables.output.tsv,
    output:
        svg = "Pipeline/SubcellularLoc/SubcellLoc.svg",
        r_is_not_a_valid_programming_language_and_no_one_should_use_this_bullshit = "Pipeline/SubcellularLoc/SubcellLocSameBullshitButCauseRCantReadDataframesItNeedsToBePreprocessed.tsv"

    script: "scripts/checkMembraneDistribution.py"



rule cpF3SubCellLoc:
    input:
        file = rules.checkProteinLocalization.output.svg,
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
        tsv = rules.joinTables.output.tsv
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
        joined = rules.joinTables.output.tsv
    output:
        trash_tsv_for_R = "Pipeline/RBP3Check/trashTSVForR.tsv",
        useful_tsv = "Pipeline/RBP3Check/rbp3joined.tsv",
        html = "Pipeline/RBP3Check/rbp3joined.html",
    script: "scripts/checkRBP3Data.py"

rule rnatypeGSEA:
    input:
        tsv = rules.joinTables.output.output_for_R_because_this_stupid_language_cannot_even_read_tsv_files
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



rule minusPVolcanoForFig1:
    input:
        tsv = rules.joinTables.output.tsv, condition="M", baseline="C",
        gsea = rules.rnatypeGSEA.output.enrichedFile, puro="minusPuromycin"
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
        tsv=rules.joinTables.output.tsv,
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



rule plotRBP3GSEA:
    input:
        df=expand(rules.rbp3GSEA.output.gsFile, puro="plusPuromycin")[0],
        info=expand(rules.rbp3GSEA.output.enrichedFile, puro="plusPuromycin")[0],
    output:
        html = "Pipeline/RBP3Check/RBP3_GSEA_plusP.html",
        svg = config["RBP3GSEA"]["file"],
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
        html = "Pipeline/Figure8/RBP3_Venn.html",
        svg = config["RBP3Venn"]["file"]
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/rbp3Venns.py"



rule joinFig8:
    input:
        files=[config[val]["file"] for key, val in config["F8"].items()]
    params:
        ges_y= 400,
        figures= config["F8"],

    output:
        svg = "Pipeline/Paper/Figure8/Figure8.svg"
    conda: "envs/svgutils.yml"
    script: "scripts/paperFigure.py"

rule Figure8toPNG:
    input:
        svg = rules.joinFig8.output.svg
    output:
        png =  "Pipeline/Paper/Figure8/Figure8.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"


rule plotReplicateCorrelations:
    input:
        tsv1 = "Pipeline/RUNS/{RUN}/PipelineData/IntermediateData/CorrData.tsv"
    params:
        layout = lambda wildcards: config[config["CorrWCMapping"][wildcards.RUN]]["layout"]
    output:
        svg = "Pipeline/Paper/Supplementary/Correlation/Run{RUN}Figure.svg"
    conda: "envs/matplotlibPlotly.yml"
    script: "scripts/plotReplicateCorrelations.py"

rule joinSupplementCorrelation:
    input:
        files=[config[val]["file"] for key, val in config["SupplementCorr"].items()]
    params:
        ges_y= 400,
        figures= config["SupplementCorr"],
    output:
        svg = "Pipeline/Paper/Supplementary/Correlation/FigureCorrelation.svg"
    conda: "envs/svgutils.yml"
    script: "scripts/paperFigure.py"


rule SupplementCorrelationtoPNG:
    input:
        svg = rules.joinSupplementCorrelation.output.svg
    output:
        png =  "Pipeline/Paper/Supplementary/Correlation/FigureCorrelation.png"
    shell: "inkscape --export-type=png --export-filename={output.png} {input.svg} --export-dpi=1200"


rule all:
    default_target: True
    input:
        tu = expand(rules.createTUTable.output, condition=["M"], baseline=["C"]),
        tu2 = expand(rules.createTUTable2.output, condition=["M"], baseline=["C"]),
        #tu3 = expand(rules.plotHistogramm.output,zip, condition=["M"], baseline=["C"]),
        stayedGO = expand(rules.GOTermEnrichmentStayed.output, condition=["M"], baseline=["C"]),
        transintramembrane = expand(rules.joinTables.output, condition=["M"], baseline=["C"]),
        size_factors = rules.third_runDESeq.output,
        #protloc = expand(rules.correlationRNAvsProteins.output, condition="M", baseline="C"),
        summary_fig = expand(rules.compareResults.output, condition="M", baseline="C"),
        summary_fig2 = expand(rules.compareTotalCell.output, condition="plusP", baseline="minusP"),
        #overlap= expand(rules.extractOverlap.output, zip, condition=config["primary"]["conditions"],baseline=config["primary"]["baselines"]),

        test = expand(rules.plotGSEA.output, condition="M", baseline="C", run_dir=[os.path.basename(config["RUNS"][run]["RUN_DIR"]) for run in config["RUNS"] if "TC" not in run]),
        subcell = rules.download_subcell.output,
        rbp3 = expand(rules.plotRBP3GSEA.output, puro=["plusPuromycin", ]),
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
        transmemgsea = expand(rules.TransmembraneGSEA.output.enriched, puro=["minusPuromycin", "plusPuromycin"], condition="M", baseline="C"),
        f5perc = config["F5Percent"]["file"],
        f5fish = config["F5FISH"]["file"],
        rbp3venn= rules.Figure8toPNG.output,
        supplement = rules.SupplementCorrelationtoPNG.output
#tpmemabundace = rules.checkTranscriptAbundanceGO.output


