#! /usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

height_full <- read.table(args[1], header=TRUE, fill=TRUE)
height_full <- height_full[is.na(height_full$f) == FALSE,]

sample_size <- floor(0.5*nrow(height_full))
set.seed(1998)

picked <- sample(seq_len(nrow(height_full)), size=sample_size)
gwas_A <- height_full[picked,]
gwas_B <- height_full[-picked,]

write.table(gwas_A,args[2],quote=FALSE,row.names=FALSE,col.names=c('id','id','height'))
write.table(gwas_B,args[3],quote=FALSE,row.names=FALSE,col.names=c('id','id','height'))
