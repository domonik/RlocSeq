

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

customGSEA <- function(table, geneList, terms, enrichedFile, gsFile, pval_cutoff) {

  TERM2GENE <- data.frame(TERM = character(0), GENE = character(0), stringsAsFactors = FALSE)
  for (term in terms) {
    mask <- table[[term]]
    print(head(mask))
    print(term)
    print(sum(mask))
    print(head(table[term]))
    term_names <- rownames(table)[mask]
    if (length(term_names) > 0) {
      TERM2GENE <- rbind(TERM2GENE, data.frame(TERM = term, GENE = term_names, stringsAsFactors = FALSE))
    }

  }
  print(head(geneList))
  print(head(TERM2GENE))

  ego <- GSEA(geneList = geneList,
              minGSSize = 1,
              maxGSSize = 10000,
              pvalueCutoff = pval_cutoff,
              verbose = FALSE,
              TERM2GENE = TERM2GENE
  )
  summary <- data.frame(ego)
  print(summary)


  if (dim(summary)[1] == 0) {
    df <- data.frame(matrix(ncol = 12, nrow = 0))
    x <- c("ONTOLOGY", "ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
    colnames(df) <- x
    summary <- df
  }

  write.table(summary, file = enrichedFile, row.names = FALSE, sep = "\t")

  m <- nrow(summary)
  if (m > 0) {
    geneSetID <- 1:m
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = ego))

  } else {
    gsdata <- data.frame(matrix(ncol = 8, nrow = 0))
    x <- c("x",	"runningScore",	"position",	"ymin",	"ymax",	"geneList",	"gene",	"Description")
    colnames(gsdata) <- x

  }

  write.table(gsdata, file = gsFile, row.names = FALSE, sep = "\t")
}



