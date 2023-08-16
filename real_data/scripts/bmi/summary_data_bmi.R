#! usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

## function to create subsetted version of results for each chromosome
sub <- function(x){
  x <- x[,c(1,2,3,9,10)]
  names(x)[1] <- 'chr'
  names(x)[2] <- 'pos'
  names(x)[3] <- 'rsid'
  names(x)[4] <- 'beta'
  names(x)[5] <- 'se'
  return(x)
}


bmiA_1 <- read.table(args[1],fill=TRUE)
bmiB_1 <- read.table(args[2],fill=TRUE)
bmiA_1 <- sub(bmiA_1)
bmiB_1 <- sub(bmiB_1)

bmiA_2 <- read.table(args[3],fill=TRUE)
bmiB_2 <- read.table(args[4],fill=TRUE)
bmiA_2 <- sub(bmiA_2)
bmiB_2 <- sub(bmiB_2)

bmiA_3 <- read.table(args[5],fill=TRUE)
bmiB_3 <- read.table(args[6],fill=TRUE)
bmiA_3 <- sub(bmiA_3)
bmiB_3 <- sub(bmiB_3)

bmiA_4 <- read.table(args[7],fill=TRUE)
bmiB_4 <- read.table(args[8],fill=TRUE)
bmiA_4 <- sub(bmiA_4)
bmiB_4 <- sub(bmiB_4)

bmiA_5 <- read.table(args[9],fill=TRUE)
bmiB_5 <- read.table(args[10],fill=TRUE)
bmiA_5 <- sub(bmiA_5)
bmiB_5 <- sub(bmiB_5)

bmiA_6 <- read.table(args[11],fill=TRUE)
bmiB_6 <- read.table(args[12],fill=TRUE)
bmiA_6 <- sub(bmiA_6)
bmiB_6 <- sub(bmiB_6)

bmiA_7 <- read.table(args[13],fill=TRUE)
bmiB_7 <- read.table(args[14],fill=TRUE)
bmiA_7 <- sub(bmiA_7)
bmiB_7 <- sub(bmiB_7)

bmiA_8 <- read.table(args[15],fill=TRUE)
bmiB_8 <- read.table(args[16],fill=TRUE)
bmiA_8 <- sub(bmiA_8)
bmiB_8 <- sub(bmiB_8)

bmiA_9 <- read.table(args[17],fill=TRUE)
bmiB_9 <- read.table(args[18],fill=TRUE)
bmiA_9 <- sub(bmiA_9)
bmiB_9 <- sub(bmiB_9)

bmiA_10 <- read.table(args[19],fill=TRUE)
bmiB_10 <- read.table(args[20],fill=TRUE)
bmiA_10 <- sub(bmiA_10)
bmiB_10 <- sub(bmiB_10)

bmiA_11 <- read.table(args[21],fill=TRUE)
bmiB_11 <- read.table(args[22],fill=TRUE)
bmiA_11 <- sub(bmiA_11)
bmiB_11 <- sub(bmiB_11)

bmiA_12 <- read.table(args[23],fill=TRUE)
bmiB_12 <- read.table(args[24],fill=TRUE)
bmiA_12 <- sub(bmiA_12)
bmiB_12 <- sub(bmiB_12)

bmiA_13 <- read.table(args[25],fill=TRUE)
bmiB_13 <- read.table(args[26],fill=TRUE)
bmiA_13 <- sub(bmiA_13)
bmiB_13 <- sub(bmiB_13)

bmiA_14 <- read.table(args[27],fill=TRUE)
bmiB_14 <- read.table(args[28],fill=TRUE)
bmiA_14 <- sub(bmiA_14)
bmiB_14 <- sub(bmiB_14)

bmiA_15 <- read.table(args[29],fill=TRUE)
bmiB_15 <- read.table(args[30],fill=TRUE)
bmiA_15 <- sub(bmiA_15)
bmiB_15 <- sub(bmiB_15)

bmiA_16 <- read.table(args[31],fill=TRUE)
bmiB_16 <- read.table(args[32],fill=TRUE)
bmiA_16 <- sub(bmiA_16)
bmiB_16 <- sub(bmiB_16)

bmiA_17 <- read.table(args[33],fill=TRUE)
bmiB_17 <- read.table(args[34],fill=TRUE)
bmiA_17 <- sub(bmiA_17)
bmiB_17 <- sub(bmiB_17)

bmiA_18 <- read.table(args[35],fill=TRUE)
bmiB_18 <- read.table(args[36],fill=TRUE)
bmiA_18 <- sub(bmiA_18)
bmiB_18 <- sub(bmiB_18)

bmiA_19 <- read.table(args[37],fill=TRUE)
bmiB_19 <- read.table(args[38],fill=TRUE)
bmiA_19 <- sub(bmiA_19)
bmiB_19 <- sub(bmiB_19)

bmiA_20 <- read.table(args[39],fill=TRUE)
bmiB_20 <- read.table(args[40],fill=TRUE)
bmiA_20 <- sub(bmiA_20)
bmiB_20 <- sub(bmiB_20)

bmiA_21 <- read.table(args[41],fill=TRUE)
bmiB_21 <- read.table(args[42],fill=TRUE)
bmiA_21 <- sub(bmiA_21)
bmiB_21 <- sub(bmiB_21)

bmiA_22 <- read.table(args[43],fill=TRUE)
bmiB_22 <- read.table(args[44],fill=TRUE)
bmiA_22 <- sub(bmiA_22)
bmiB_22 <- sub(bmiB_22)


bmiA <- rbind(bmiA_1,bmiA_2,bmiA_3,bmiA_4,bmiA_5,bmiA_6,bmiA_7,bmiA_8,bmiA_9,bmiA_10,bmiA_11,bmiA_12,bmiA_13,bmiA_14,bmiA_15,bmiA_16,bmiA_17,bmiA_18,bmiA_19,bmiA_20,bmiA_21,bmiA_22)
bmiB <- rbind(bmiB_1,bmiB_2,bmiB_3,bmiB_4,bmiB_5,bmiB_6,bmiB_7,bmiB_8,bmiB_9,bmiB_10,bmiB_11,bmiB_12,bmiB_13,bmiB_14,bmiB_15,bmiB_16,bmiB_17,bmiB_18,bmiB_19,bmiB_20,bmiB_21,bmiB_22)


if(sum(bmiA$rsid != bmiB$rsid)==0){bmiA$beta_rep <- bmiB$beta}

write.table(bmiA, args[45], quote=FALSE, row.names=FALSE)
