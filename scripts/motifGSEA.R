source("scripts/customGSEA.R")

detable <- snakemake@input[["tsv"]]
detable <- read.table(detable, sep="\t", header=TRUE)
rownames(detable) <- detable$locus_tag
detable <- detable[!is.na(detable$log2FoldChange_plusPuromycin), ]

df_sorted <- detable[order(detable$log2FoldChange_plusPuromycin, decreasing = TRUE), ]
geneList <- df_sorted$log2FoldChange_plusPuromycin
names(geneList) <- rownames(df_sorted)
terms <- grep("Motif", colnames(df_sorted), value = TRUE)

# Convert these columns to logical (boolean)
df_sorted[, terms] <- lapply(df_sorted[, terms], function(x) x == "True") # why the fuck is this necessary

efile <- snakemake@output[["enrichedFile"]]
gsfile <- snakemake@output[["gsFile"]]

customGSEA(
  df_sorted,
  geneList = geneList,
  terms = terms,
  enrichedFile = efile,
  gsFile = gsfile,
  pval_cutoff = 1

)


