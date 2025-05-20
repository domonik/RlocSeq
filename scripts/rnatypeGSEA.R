source("scripts/customGSEA.R")

detable <- snakemake@input[["tsv"]]
detable <- read.table(detable, sep="\t", header=TRUE, check.names = FALSE)
print(head(detable))

rownames(detable) <- detable$locus_tag
fc <- paste0("log2FoldChange_", snakemake@wildcards[["puro"]])
print(fc)
detable <- detable[!is.na(detable[[fc]]), ]
df_sorted <- detable[order(detable[[fc]], decreasing = TRUE), ]

geneList <- df_sorted[[fc]]
names(geneList) <- rownames(df_sorted)
terms <- snakemake@params[["terms"]]
print(terms)
print(typeof(terms))
print(head(df_sorted))
terms <- intersect(terms, colnames(df_sorted))
# Convert these columns to logical (boolean)
df_sorted[, terms] <- lapply(df_sorted[, terms], function(x) x == "True") # why the fuck is this necessary

# Convert these columns to logical (boolean)
# why the fuck is this necessary
efile <- snakemake@output[["enrichedFile"]]
gsfile <- snakemake@output[["gsFile"]]

customGSEA(
  df_sorted,
  geneList = geneList,
  terms = terms,
  enrichedFile = efile,
  gsFile = gsfile,
  pval_cutoff = 0.05,
)


