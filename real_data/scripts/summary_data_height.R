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


heightA_1 <- read.table(args[1],fill=TRUE)
heightB_1 <- read.table(args[2],fill=TRUE)
heightA_1 <- sub(heightA_1)
heightB_1 <- sub(heightB_1)

heightA_2 <- read.table(args[3],fill=TRUE)
heightB_2 <- read.table(args[4],fill=TRUE)
heightA_2 <- sub(heightA_2)
heightB_2 <- sub(heightB_2)

heightA_3 <- read.table(args[5],fill=TRUE)
heightB_3 <- read.table(args[6],fill=TRUE)
heightA_3 <- sub(heightA_3)
heightB_3 <- sub(heightB_3)

heightA_4 <- read.table(args[7],fill=TRUE)
heightB_4 <- read.table(args[8],fill=TRUE)
heightA_4 <- sub(heightA_4)
heightB_4 <- sub(heightB_4)

heightA_5 <- read.table(args[9],fill=TRUE)
heightB_5 <- read.table(args[10],fill=TRUE)
heightA_5 <- sub(heightA_5)
heightB_5 <- sub(heightB_5)

heightA_6 <- read.table(args[11],fill=TRUE)
heightB_6 <- read.table(args[12],fill=TRUE)
heightA_6 <- sub(heightA_6)
heightB_6 <- sub(heightB_6)

heightA_7 <- read.table(args[13],fill=TRUE)
heightB_7 <- read.table(args[14],fill=TRUE)
heightA_7 <- sub(heightA_7)
heightB_7 <- sub(heightB_7)

heightA_8 <- read.table(args[15],fill=TRUE)
heightB_8 <- read.table(args[16],fill=TRUE)
heightA_8 <- sub(heightA_8)
heightB_8 <- sub(heightB_8)

heightA_9 <- read.table(args[17],fill=TRUE)
heightB_9 <- read.table(args[18],fill=TRUE)
heightA_9 <- sub(heightA_9)
heightB_9 <- sub(heightB_9)

heightA_10 <- read.table(args[19],fill=TRUE)
heightB_10 <- read.table(args[20],fill=TRUE)
heightA_10 <- sub(heightA_10)
heightB_10 <- sub(heightB_10)

heightA_11 <- read.table(args[21],fill=TRUE)
heightB_11 <- read.table(args[22],fill=TRUE)
heightA_11 <- sub(heightA_11)
heightB_11 <- sub(heightB_11)

heightA_12 <- read.table(args[23],fill=TRUE)
heightB_12 <- read.table(args[24],fill=TRUE)
heightA_12 <- sub(heightA_12)
heightB_12 <- sub(heightB_12)

heightA_13 <- read.table(args[25],fill=TRUE)
heightB_13 <- read.table(args[26],fill=TRUE)
heightA_13 <- sub(heightA_13)
heightB_13 <- sub(heightB_13)

heightA_14 <- read.table(args[27],fill=TRUE)
heightB_14 <- read.table(args[28],fill=TRUE)
heightA_14 <- sub(heightA_14)
heightB_14 <- sub(heightB_14)

heightA_15 <- read.table(args[29],fill=TRUE)
heightB_15 <- read.table(args[30],fill=TRUE)
heightA_15 <- sub(heightA_15)
heightB_15 <- sub(heightB_15)

heightA_16 <- read.table(args[31],fill=TRUE)
heightB_16 <- read.table(args[32],fill=TRUE)
heightA_16 <- sub(heightA_16)
heightB_16 <- sub(heightB_16)

heightA_17 <- read.table(args[33],fill=TRUE)
heightB_17 <- read.table(args[34],fill=TRUE)
heightA_17 <- sub(heightA_17)
heightB_17 <- sub(heightB_17)

heightA_18 <- read.table(args[35],fill=TRUE)
heightB_18 <- read.table(args[36],fill=TRUE)
heightA_18 <- sub(heightA_18)
heightB_18 <- sub(heightB_18)

heightA_19 <- read.table(args[37],fill=TRUE)
heightB_19 <- read.table(args[38],fill=TRUE)
heightA_19 <- sub(heightA_19)
heightB_19 <- sub(heightB_19)

heightA_20 <- read.table(args[39],fill=TRUE)
heightB_20 <- read.table(args[40],fill=TRUE)
heightA_20 <- sub(heightA_20)
heightB_20 <- sub(heightB_20)

heightA_21 <- read.table(args[41],fill=TRUE)
heightB_21 <- read.table(args[42],fill=TRUE)
heightA_21 <- sub(heightA_21)
heightB_21 <- sub(heightB_21)

heightA_22 <- read.table(args[43],fill=TRUE)
heightB_22 <- read.table(args[44],fill=TRUE)
heightA_22 <- sub(heightA_22)
heightB_22 <- sub(heightB_22)


heightA <- rbind(heightA_1,heightA_2,heightA_3,heightA_4,heightA_5,heightA_6,heightA_7,heightA_8,heightA_9,heightA_10,heightA_11,heightA_12,heightA_13,heightA_14,heightA_15,heightA_16,heightA_17,heightA_18,heightA_19,heightA_20,heightA_21,heightA_22)
heightB <- rbind(heightB_1,heightB_2,heightB_3,heightB_4,heightB_5,heightB_6,heightB_7,heightB_8,heightB_9,heightB_10,heightB_11,heightB_12,heightB_13,heightB_14,heightB_15,heightB_16,heightB_17,heightB_18,heightB_19,heightB_20,heightB_21,heightB_22)


if(sum(heightA$rsid != heightB$rsid)==0){heightA$beta_rep <- heightB$beta}

write.table(heightA, args[45], quote=FALSE, row.names=FALSE)
