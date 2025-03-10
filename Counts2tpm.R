## Counts to TPM

library(IOBR)
setwd("~/Documents/practiques/data")

data <- read.csv("unproc_rnaseq.csv", row.names = 1)
head(df, 10)

data_tpm <- count2tpm(countMat = data, source = "local", idType = "Entrez")
head(data_tpm)

write.csv(data_tpm, "TPM_results.csv", row.names = TRUE)
