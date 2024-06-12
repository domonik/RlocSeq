

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  require(svglite)
})
gseaScores <- getFromNamespace("gseaScores", "DOSE")

gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    } else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}



table <- snakemake@input[["deseq_results"]]
table <- read.table(table, sep="\t", header=TRUE)
table <- table[!is.na(table$log2FoldChange), ]
df_sorted <- table[order(table$log2FoldChange, decreasing = TRUE), ]

# Step 2: Identify row names containing "5UTR" or "3UTR"
rownames_table <- rownames(df_sorted)
utr5_names <- rownames_table[grep("5UTR", rownames_table)]
utr3_names <- rownames_table[grep("3UTR", rownames_table)]
nonCoding_names <- rownames_table[grep("ncr|ncl", rownames_table)]
trna_names <- rownames_table[grep("6803t", rownames_table)]
asrna_names <- rownames_table[grep("-as", rownames_table)]
crisprrna_names <- rownames_table[grep("CRISPR", rownames_table)]
# Exclude row names that contain "-as"
mrna_names <- rownames_table[grep("slr|sll|ssr|ssl", rownames_table)]
mrna_names <- mrna_names[!grepl("-as", mrna_names)]
mrna_names <- mrna_names[!grepl("UTR", mrna_names)]

random <- sample(rownames_table, 2000)


TERM2GENE <- data.frame(TERM = character(0), GENE = character(0), stringsAsFactors = FALSE)
if (length(utr5_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "5UTR", GENE = utr5_names, stringsAsFactors = FALSE))
}
if (length(nonCoding_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "ncRNA", GENE = nonCoding_names, stringsAsFactors = FALSE))
}
if (length(trna_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "tRNA", GENE = trna_names, stringsAsFactors = FALSE))
}

if (length(mrna_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "mRNA", GENE = mrna_names, stringsAsFactors = FALSE))
}

if (length(mrna_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "random", GENE = random, stringsAsFactors = FALSE))
}

# Add entries for 3UTR
if (length(utr3_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "3UTR", GENE = utr3_names, stringsAsFactors = FALSE))
}

if (length(asrna_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "asRNA", GENE = asrna_names, stringsAsFactors = FALSE))
}
if (length(crisprrna_names) > 0) {
  TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = "CRISPR", GENE = crisprrna_names, stringsAsFactors = FALSE))
}



geneList <- df_sorted$log2FoldChange
names(geneList) <- rownames(df_sorted)

ego <- GSEA(geneList = geneList,
              minGSSize = 5,
              maxGSSize = 10000,
              pvalueCutoff = 1,
              verbose = FALSE,
              TERM2GENE=TERM2GENE
)
summary <- data.frame(ego)


if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
  colnames(df) <- x
  summary <- df
}

write.table(summary , file = snakemake@output[["enriched"]], row.names=FALSE, sep="\t")

m <- nrow(summary)
geneSetID <- 1:m
gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = ego))
write.table(gsdata , file = snakemake@output[["gsdata"]], row.names=FALSE, sep="\t")

