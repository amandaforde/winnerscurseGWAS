#! usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

## function to create subsetted version of results for each chromosome
sub <- function(x){
  x <- x[,c(1,2,3,10,11)]
  names(x)[1] <- 'chr'
  names(x)[2] <- 'pos'
  names(x)[3] <- 'rsid'
  names(x)[4] <- 'beta'
  names(x)[5] <- 'se'
  return(x)
}


T2DA_1 <- read.table(args[1],fill=TRUE)
T2DB_1 <- read.table(args[2],fill=TRUE)
T2DA_1 <- sub(T2DA_1)
T2DB_1 <- sub(T2DB_1)

T2DA_2 <- read.table(args[3],fill=TRUE)
T2DB_2 <- read.table(args[4],fill=TRUE)
T2DA_2 <- sub(T2DA_2)
T2DB_2 <- sub(T2DB_2)

T2DA_3 <- read.table(args[5],fill=TRUE)
T2DB_3 <- read.table(args[6],fill=TRUE)
T2DA_3 <- sub(T2DA_3)
T2DB_3 <- sub(T2DB_3)

T2DA_4 <- read.table(args[7],fill=TRUE)
T2DB_4 <- read.table(args[8],fill=TRUE)
T2DA_4 <- sub(T2DA_4)
T2DB_4 <- sub(T2DB_4)

T2DA_5 <- read.table(args[9],fill=TRUE)
T2DB_5 <- read.table(args[10],fill=TRUE)
T2DA_5 <- sub(T2DA_5)
T2DB_5 <- sub(T2DB_5)

T2DA_6 <- read.table(args[11],fill=TRUE)
T2DB_6 <- read.table(args[12],fill=TRUE)
T2DA_6 <- sub(T2DA_6)
T2DB_6 <- sub(T2DB_6)

T2DA_7 <- read.table(args[13],fill=TRUE)
T2DB_7 <- read.table(args[14],fill=TRUE)
T2DA_7 <- sub(T2DA_7)
T2DB_7 <- sub(T2DB_7)

T2DA_8 <- read.table(args[15],fill=TRUE)
T2DB_8 <- read.table(args[16],fill=TRUE)
T2DA_8 <- sub(T2DA_8)
T2DB_8 <- sub(T2DB_8)

T2DA_9 <- read.table(args[17],fill=TRUE)
T2DB_9 <- read.table(args[18],fill=TRUE)
T2DA_9 <- sub(T2DA_9)
T2DB_9 <- sub(T2DB_9)

T2DA_10 <- read.table(args[19],fill=TRUE)
T2DB_10 <- read.table(args[20],fill=TRUE)
T2DA_10 <- sub(T2DA_10)
T2DB_10 <- sub(T2DB_10)

T2DA_11 <- read.table(args[21],fill=TRUE)
T2DB_11 <- read.table(args[22],fill=TRUE)
T2DA_11 <- sub(T2DA_11)
T2DB_11 <- sub(T2DB_11)

T2DA_12 <- read.table(args[23],fill=TRUE)
T2DB_12 <- read.table(args[24],fill=TRUE)
T2DA_12 <- sub(T2DA_12)
T2DB_12 <- sub(T2DB_12)

T2DA_13 <- read.table(args[25],fill=TRUE)
T2DB_13 <- read.table(args[26],fill=TRUE)
T2DA_13 <- sub(T2DA_13)
T2DB_13 <- sub(T2DB_13)

T2DA_14 <- read.table(args[27],fill=TRUE)
T2DB_14 <- read.table(args[28],fill=TRUE)
T2DA_14 <- sub(T2DA_14)
T2DB_14 <- sub(T2DB_14)

T2DA_15 <- read.table(args[29],fill=TRUE)
T2DB_15 <- read.table(args[30],fill=TRUE)
T2DA_15 <- sub(T2DA_15)
T2DB_15 <- sub(T2DB_15)

T2DA_16 <- read.table(args[31],fill=TRUE)
T2DB_16 <- read.table(args[32],fill=TRUE)
T2DA_16 <- sub(T2DA_16)
T2DB_16 <- sub(T2DB_16)

T2DA_17 <- read.table(args[33],fill=TRUE)
T2DB_17 <- read.table(args[34],fill=TRUE)
T2DA_17 <- sub(T2DA_17)
T2DB_17 <- sub(T2DB_17)

T2DA_18 <- read.table(args[35],fill=TRUE)
T2DB_18 <- read.table(args[36],fill=TRUE)
T2DA_18 <- sub(T2DA_18)
T2DB_18 <- sub(T2DB_18)

T2DA_19 <- read.table(args[37],fill=TRUE)
T2DB_19 <- read.table(args[38],fill=TRUE)
T2DA_19 <- sub(T2DA_19)
T2DB_19 <- sub(T2DB_19)

T2DA_20 <- read.table(args[39],fill=TRUE)
T2DB_20 <- read.table(args[40],fill=TRUE)
T2DA_20 <- sub(T2DA_20)
T2DB_20 <- sub(T2DB_20)

T2DA_21 <- read.table(args[41],fill=TRUE)
T2DB_21 <- read.table(args[42],fill=TRUE)
T2DA_21 <- sub(T2DA_21)
T2DB_21 <- sub(T2DB_21)

T2DA_22 <- read.table(args[43],fill=TRUE)
T2DB_22 <- read.table(args[44],fill=TRUE)
T2DA_22 <- sub(T2DA_22)
T2DB_22 <- sub(T2DB_22)


T2DA <- rbind(T2DA_1,T2DA_2,T2DA_3,T2DA_4,T2DA_5,T2DA_6,T2DA_7,T2DA_8,T2DA_9,T2DA_10,T2DA_11,T2DA_12,T2DA_13,T2DA_14,T2DA_15,T2DA_16,T2DA_17,T2DA_18,T2DA_19,T2DA_20,T2DA_21,T2DA_22)
T2DB <- rbind(T2DB_1,T2DB_2,T2DB_3,T2DB_4,T2DB_5,T2DB_6,T2DB_7,T2DB_8,T2DB_9,T2DB_10,T2DB_11,T2DB_12,T2DB_13,T2DB_14,T2DB_15,T2DB_16,T2DB_17,T2DB_18,T2DB_19,T2DB_20,T2DB_21,T2DB_22)


if(sum(T2DA$rsid != T2DB$rsid)==0){T2DA$beta_rep <- T2DB$beta}

write.table(T2DA, args[45], quote=FALSE, row.names=FALSE)
