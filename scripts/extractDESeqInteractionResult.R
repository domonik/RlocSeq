suppressPackageStartupMessages({
  require(DESeq2)
})

file <- snakemake@input[["result"]]
contrast_coefficients <- unlist(
  snakemake@params[["contrast_coefficients"]],
  use.names = FALSE
)

load(file = file)

available_coefficients <- resultsNames(dds)
missing_coefficients <- setdiff(contrast_coefficients, available_coefficients)
if (length(missing_coefficients) > 0) {
  stop(
    "Could not construct the interaction contrast. Missing DESeq2 coefficient(s): ",
    paste(missing_coefficients, collapse = ", "),
    ". Available coefficients: ",
    paste(available_coefficients, collapse = ", ")
  )
}

# The interaction tests whether membrane-versus-cytosol enrichment changes
# under puromycin, cancelling a treatment abundance effect shared by both
# fractions. Biological samples are nested within treatment in the fitted
# design, so equally numbered replicates are not paired across treatments.
# SynDESeq's normal-prior shrinkage is omitted because DESeq2 does not implement
# that shrinkage for this design.
if (length(contrast_coefficients) == 1) {
  result_table <- results(dds, name = contrast_coefficients, alpha = 0.05)
} else {
  result_table <- results(
    dds,
    contrast = list(contrast_coefficients),
    alpha = 0.05
  )
}

write.table(result_table, snakemake@output[["result_table"]], sep = "\t")
