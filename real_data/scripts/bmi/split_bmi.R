#! /usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

bmi_full <- read.table(args[1], header=TRUE, fill=TRUE)
bmi_full <- bmi_full[is.na(bmi_full$f) == FALSE,]

sample_size <- floor(0.5*nrow(bmi_full))
set.seed(1998)

picked <- sample(seq_len(nrow(bmi_full)), size=sample_size)
gwas_A <- bmi_full[picked,]
gwas_B <- bmi_full[-picked,]

write.table(gwas_A,args[2],quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))
write.table(gwas_B,args[3],quote=FALSE,row.names=FALSE,col.names=c('id','id','bmi'))



