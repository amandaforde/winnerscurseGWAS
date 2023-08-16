#! /usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

T2D_full <- read.table(args[1], header=TRUE, fill=TRUE)
T2D_full <- T2D_full[is.na(T2D_full$T2D) == FALSE,]

samples_size <- floor(0.5*nrow(T2D_full))
set.seed(1998)

picked <- sample(seq_len(nrow(T2D_full)), size=sample_size)
gwas_A <- T2D_full[picked,]
gwas_B <- T2D_full[-picked,]

gwas_A$T2D <- gwas_A$T2D + 1
gwas_B$T2D <- gwas_B$T2D + 1

write.table(gwas_A,args[2],quote=FALSE,row.names=FALSE,col.names=c('id','id','T2D'))
write.table(gwas_B,args[3],quote=FALSE,row.names=FALSE,col.names=c('id','id','T2D'))


