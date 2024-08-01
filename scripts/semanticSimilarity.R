library(Rcpp)
library(GOSemSim)
package <- list.files(snakemake@input[["annotation_db"]])[1]
library(basename(package), character.only = TRUE, lib.loc=snakemake@config[["Rlib"]])


gosubcat <- snakemake@wildcards[["subcat"]]

file <- snakemake@input[["enrichData"]]
table <- read.table(file ,sep="\t",header=TRUE)

subtable <- table[table$ONTOLOGY == gosubcat,]
semdata <- godata(basename(package), ont=gosubcat, computeIC=T, keytype="GID")
semsim <- mgoSim(subtable$ID, subtable$ID, semData=semdata, measure="Rel", combine=NULL)
semsim <- as.data.frame(semsim)


write.table(semsim , file = snakemake@output[["table"]], row.names=T, sep="\t")
