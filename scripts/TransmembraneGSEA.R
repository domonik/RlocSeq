source("scripts/customGSEA.R")

detable <- snakemake@input[["tsv"]]
detable <- read.table(detable, sep="\t", header=TRUE, check.names = FALSE)
print(head(detable))

rownames(detable) <- detable$locus_tag
fc <- paste0("log2FoldChange_", snakemake@wildcards[["puro"]])
padj <- paste0("padj_", snakemake@wildcards[["puro"]])
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
print(head(df_sorted[["Transmembrane"]], n=50))

df_sorted[, terms] <- lapply(df_sorted[, terms, drop = FALSE], function(x) {
  as.logical(trimws(as.character(x)))
})
print(head(df_sorted[["Transmembrane"]], n=50))
df_sorted$NonTransmembrane <- !df_sorted$Transmembrane
efile <- snakemake@output[["enriched"]]
terms <- snakemake@params[["terms"]]


TERM2GENE <- data.frame(TERM = character(0), GENE = character(0), stringsAsFactors = FALSE)
for (term in terms) {
  print("term")
  print(term)
  mask <- df_sorted[[term]]

  term_names <- rownames(df_sorted)[mask]
  if (length(term_names) > 0) {
    TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = term, GENE = term_names, stringsAsFactors = FALSE))
  }

}

down <- df_sorted[[fc]] <= -0.8 & df_sorted[[padj]] < 0.05
down <- df_sorted[down, ]
print(TERM2GENE)
print(head(down))
print((dim(down)))
summary <- enricher(
  gene          = as.character(down$locus_tag),
  TERM2GENE     = TERM2GENE, # data.frame with TERM and GENE
  universe      = as.character(df_sorted$locus_tag),
  maxGSSize = 2000,
                    pvalueCutoff = 1,
                    qvalueCutoff =1
)
summary <- data.frame(summary)
write.table(summary, file = efile, row.names = FALSE, sep = "\t")




