## WINNER'S CURSE REAL DATA FIGURES:

## This script provides the code used to produce the figures in the manuscript
## concerning real data.

## Outputs:
## 1. Fig2.eps  
## 2. Fig3.eps 
## 3. S16_Fig.tiff - S22_Fig.tiff

## Load required packages:
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(qqman)
library("RColorBrewer")
col <- brewer.pal(8,"Dark2")
col1 <- brewer.pal(11,"RdYlBu")

## Load data sets:
summary_data_bmi_1 <-  read.table('real_data/data/summary_data_bmi_1.txt',header=TRUE)
summary_data_bmi_2 <-  read.table('real_data/data/summary_data_bmi_2.txt',header=TRUE)

summary_data_T2D_1 <-  read.table('real_data/data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)
summary_data_T2D_2 <-  read.table('real_data/data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)

summary_data_height_1 <-  read.table('real_data/data/summary_data_height_1.txt',header=TRUE)
summary_data_height_2 <-  read.table('real_data/data/summary_data_height_2.txt',header=TRUE)

################################################################################

## Fig 2. Estimated MSE of significant SNPs at threshold 5 × 10-8  for each
## method and data set.

mse_5e_8 <- read.csv("real_data/results/mse_5e_8.txt")
new_mse_5e_8 <- pivot_longer(mse_5e_8, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_8$Method <- factor(new_mse_5e_8$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_8$GWAS <- factor(new_mse_5e_8$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_1 <- ggplot(new_mse_5e_8,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-8))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_8[new_mse_5e_8$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_1
## Save as 'figures/Fig2.eps' with dimensions 800 x 1100 pixels.

################################################################################

## Fig 3. Estimated MSE of significant SNPs at threshold 5 × 10-4 for each
## method and data set.

mse_5e_4 <- read.csv("real_data/results/mse_5e_4.txt")
new_mse_5e_4 <- pivot_longer(mse_5e_4, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_4$Method <- factor(new_mse_5e_4$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_4$GWAS <- factor(new_mse_5e_4$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_2 <- ggplot(new_mse_5e_4,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-4))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_4[new_mse_5e_4$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_2
## Save as 'figures/Fig3.eps' with dimensions 800 x 1100 pixels.

################################################################################

## S16 Fig. Z-statistics plotted against bias for each real data set.

# BMI GWAS 1
out <- summary_data_bmi_1
title <- "BMI GWAS 1"
title_col <- col[1]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
BMI1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                                      aes(x=z,y=bias),
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )

# BMI GWAS 2
out <- summary_data_bmi_2
title <- "BMI GWAS 2"
title_col <- col[2]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
BMI2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                                      aes(x=z,y=bias),
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )

# T2D GWAS 1
out <- summary_data_T2D_1
title <- "T2D GWAS 1"
title_col <- col[3]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
T2D1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                                      aes(x=z,y=bias),
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )

# T2D GWAS 2
out <- summary_data_T2D_2
title <- "T2D GWAS 2"
title_col <- col[4]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
T2D2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                                      aes(x=z,y=bias),
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )


# Height GWAS 1
out <- summary_data_height_1
title <- "Height GWAS 1"
title_col <- col[5]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
Height1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                                         aes(x=z,y=bias),
                                                                         color='grey40',
                                                                         size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )


# Height GWAS 2
out <- summary_data_height_2
title <- "Height GWAS 2"
title_col <- col[6]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
Height2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout,
                                                                         aes(x=z,y=bias),
                                                                         color='grey40',
                                                                         size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )

# Combine above six plots into a figure together
plot_3 <- ggarrange(BMI1, BMI2, T2D1, T2D2, Height1, Height2, labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

plot_3
## Save as 'figures/S16_Fig.tiff' with dimensions 1200 x 1200 pixels.

################################################################################

## S17 Fig. Average bias of significant SNPs with positive association estimates
## at threshold 5 × 10-8 for each method and data set.

bias_5e_8_up <- read.csv("real_data/results/bias_5e_8_positive.txt")
new_bias_5e_8_up <- pivot_longer(bias_5e_8_up, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_8_up$Method <- factor(new_bias_5e_8_up$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_8_up$GWAS <- factor(new_bias_5e_8_up$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_4 <- ggplot(new_bias_5e_8_up,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-8, " for positive SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_bias_5e_8_up[new_bias_5e_8_up$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_4
## Save as 'figures/S17_Fig.tiff' with dimensions 800 x 1100 pixels.

################################################################################

## S18 Fig. Average bias of significant SNPs with negative association estimates
## at threshold 5 × 10-8 for each method and data set.

bias_5e_8_down <- read.csv("real_data/results/bias_5e_8_negative.txt")
new_bias_5e_8_down <- pivot_longer(bias_5e_8_down, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_8_down$Method <- factor(new_bias_5e_8_down$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_8_down$GWAS <- factor(new_bias_5e_8_down$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_5 <- ggplot(new_bias_5e_8_down,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-8, " for negative SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_bias_5e_8_down[new_bias_5e_8_down$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_5
## Save as 'figures/S18_Fig.tiff' with dimensions 800 x 1100 pixels.

################################################################################

## S19 Fig. Average bias of significant SNPs with positive association estimates
## at threshold 5 × 10-4 for each method and data set.

bias_5e_4_up <- read.csv("real_data/results/bias_5e_4_positive.txt")
new_bias_5e_4_up <- pivot_longer(bias_5e_4_up, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_4_up$Method <- factor(new_bias_5e_4_up$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_4_up$GWAS <- factor(new_bias_5e_4_up$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_6 <- ggplot(new_bias_5e_4_up,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-4, " for positive SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_bias_5e_4_up[new_bias_5e_4_up$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_6
## Save as 'figures/S19_Fig.tiff' with dimensions 800 x 1100 pixels.

################################################################################

## S20 Fig. Average bias of significant SNPs with negative association estimates
## at threshold 5 × 10-4 for each method and data set.

bias_5e_4_down <- read.csv("real_data/results/bias_5e_4_negative.txt")
new_bias_5e_4_down <- pivot_longer(bias_5e_4_down, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_4_down$Method <- factor(new_bias_5e_4_down$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_4_down$GWAS <- factor(new_bias_5e_4_down$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_7 <- ggplot(new_bias_5e_4_down,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-4, " for negative SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_bias_5e_4_down[new_bias_5e_4_down$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_7
## Save as 'figures/S20_Fig.tiff' with dimensions 800 x 1100 pixels.

################################################################################

## S21 Fig. Estimated MSE of significant SNPs at threshold 5 × 10-8 for each
## method and pruned data set.

mse_5e_8 <- read.csv("real_data/results/pruned_mse_5e_8.txt")
new_mse_5e_8 <- pivot_longer(mse_5e_8, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_8$Method <- factor(new_mse_5e_8$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_8$GWAS <- factor(new_mse_5e_8$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_8 <- ggplot(new_mse_5e_8,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-8))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_8[new_mse_5e_8$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_8
## Save as 'figures/S21_Fig.tiff' with dimensions 800 x 1100 pixels.

################################################################################

## S22 Fig. Estimated MSE of significant SNPs at threshold 5 × 10-4 for each
## method and pruned data set.

mse_5e_4 <- read.csv("real_data/results/pruned_mse_5e_4.txt")
new_mse_5e_4 <- pivot_longer(mse_5e_4, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_4$Method <- factor(new_mse_5e_4$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_4$GWAS <- factor(new_mse_5e_4$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))

plot_9 <- ggplot(new_mse_5e_4,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-4))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_4[new_mse_5e_4$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_9
## Save as 'figures/S22_Fig.tiff' with dimensions 800 x 1100 pixels.

################################################################################











