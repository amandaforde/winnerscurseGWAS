#! usr/bin/env Rscript

args=commandArgs(trailinOnly=TRUE)

pruned <- read.table(args[1], header=FALSE)

for (i in 2:22){
 pruned_new <- read.table(args[i],header=FALSE)
 pruned <- rbind(pruned, pruned_new)
}

write.table(pruned,args[23],quote=FALSE,row.names=FALSE)


