# Review and further developments in statistical corrections for Winner’s Curse in genetic association studies

This repository contains the scripts used to provide results and figures for the manuscript, ***"Review and further developments in statistical corrections for Winner’s Curse in genetic association studies"***. In our work, we review existing Winner’s Curse correction methods which require only GWAS summary statistics in order to make adjustments. In addition, we have proposed modifications to improve existing methods and proposed a novel approach which uses the parametric bootstrap. We evaluate and compare these methods, first using a **simulation study** with a wide variety of simulated data sets and then, by means of an **empirical analysis** using real data sets for three different traits. 

Scripts related to simulations are contained in `simulations/scripts` while scripts used for the empirical analysis are contained in `real_data/scripts`. 

&nbsp;

## Simulation study

**sims_pipeline.R:** This script allows us to run the entire simulation study which evaluates and compares various winner's curse correction methods, from beginning to end. ***Note that seeds have been included in each script separately.*** 

**useful_funs.R:** This script provides all functions required for simulations.

**nsig_prop_bias_LD.R:** This script performs a preliminary investigation of the simulated sets of summary statistics, in which a simple *correlation structure* has been imposed on independent blocks of 100 SNPs. The simulated data sets also correspond to a quantitative trait with a normal effect size distribution. At two significance thresholds, 5 <span>&#215;</span> 10<sup>-8</sup> and 5 <span>&#215;</span> 10<sup>-4</sup>, this script obtains the number of significant SNPs, the proportion of these SNPs for which their association estimate is more extreme than their true effect size, the proportion of these SNPs which are significantly overestimated and the mean square error (MSE) of significant SNPs for each simulated data set.

- ***Output:*** nsig_prop_bias_5e_8_LD_all.csv, nsig_prop_bias_5e-8_LD.csv, nsig_prop_bias_5e_4_LD_all.csv, nsig_prop_bias_5e-4_LD.csv


&nbsp;

### Evaluating methods 

**sims_LD.R:** This script evaluates and compares several winner's curse correction methods using simulated sets of summary statistics, in which a simple *correlation structure* has been imposed on independent blocks of 100 SNPs. The simulated data sets also correspond to a quantitative trait with a normal effect size distribution. At two significance thresholds, 5 <span>&#215;</span> 10<sup>-8</sup> and 5 <span>&#215;</span> 10<sup>-4</sup>, this script compares the various methods using different evaluation metrics.

- ***Output:*** norm_5e-8_100sim_LD_ave.csv, norm_5e-8_100sim_LD_all.csv, norm_5e-4_100sim_LD_ave.csv, norm_5e-4_100sim_LD_all.csv

**sims_ind_1.R:** This script evaluates and compares several winner's curse correction methods using simulated sets of summary statistics, in which SNPs are *independent*.The simulated data sets also correspond to a quantitative trait with a normal effect size distribution. At two significance thresholds, 5 <span>&#215;</span> 10<sup>-8</sup> and 5 <span>&#215;</span> 10<sup>-4</sup>, this script compares the various methods using different evaluation metrics.

- ***Output:*** norm_5e-8_100sim_ave.csv, norm_5e-8_100sim_all.csv, norm_5e-4_100sim_ave.csv, norm_5e-4_100sim_all.csv

**sims_ind_2.R:** This script evaluates and compares several winner's curse correction methods using simulated sets of summary statistics, in which SNPs are *independent*. The simulated data sets correspond to three different settings: a quantitative trait with bimodal effect size distribution, a quantitative trait with skewed effect size distribution, and a binary trait with normal effect size distribution. 5 <span>&#215;</span> 10<sup>-8</sup> and 5 <span>&#215;</span> 10<sup>-4</sup>, this script compares the various methods using different evaluation metrics.

- ***Output:*** bim_5e-8_100sim_ave.csv, bim_5e-8_100sim_all.csv, bim_5e-4_100sim_ave.csv, bim_5e-4_100sim_all.csv, skew_5e-8_100sim_ave.csv, skew_5e-8_100sim_all.csv, skew_5e-4_100sim_ave.csv, skew_5e-4_100sim_all.csv, bin_5e-8_100sim_ave.csv, bin_5e-8_100sim_all.csv, bin_5e-4_100sim_ave.csv, bin_5e-4_100sim_all.csv


&nbsp;

## Empirical analysis

### Quality control

The genotypic data used was collected, processed and imputed by UK Biobank (UKBB) ([http://www.ukbiobank.ac.uk/](http://www.ukbiobank.ac.uk/)). The required quality control steps were implemented using the code available at [https://github.com/coggene/UK-Biobank-QC](https://github.com/coggene/UK-Biobank-QC).

These steps included the removal of both related individuals and those that were of non-European ancestry. These non-European samples were identified by principal component analysis (PCA) using 1000 Genomes Project (1KGP) data. Furthermore, samples which had been identified as outliers with respect to heterozygosity and missingness, together with samples with discordant sex information and those suffering from chromosomal aneuploidy, were also discarded. The total number of samples remaining after the execution of these quality control steps were 332,618, 333,642 and 332,927 for BMI, T2D and height, respectively. With respect to variants, only those with an information score greater than 0.8, a minor allele frequency greater than 0.01, a genotyping rate of at least 98% and those that passed the Hardy-Weinberg test at the specified significance threshold of 1 <span>&#215;</span> 10<sup>-8</sup> were included. This resulted in a total of 7,915,560 SNPs that were considered suitable for our analyses. 


&nbsp;

### Generation of summary statistics

In the following, we detail the procedure followed in order to obtain two sets of GWAS summary statistics for each trait. We provide a detailed description with respect to one trait, namely BMI. The two sets of summary statistics for height and T2D were obtained in an identical manner using their respective scripts. The scripts mentioned below can be found in `real_data/scripts/bmi`.

Note that the QC process provided us with the following files: bmi4plink.txt, T2D4plink.txt, height4plink.txt, covariates.txt and \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome.

The following code was run:

`Rscript split_bmi.R bmi4plink.txt bmi_gwasA.txt bmi_gwasB.txt`   

`mkdir {bmi_gwasA,bmi_gwasB,results_bmi}`   

`sbatch parallelA_bmi.sh`   

`sbatch parallelB_bmi.sh`   

`sbatch summary_data_bmi_1.sh`   

`sbatch summary_data_bmi_2.sh` 

&nbsp;

#### Step 1: Splitting the data set  


**split_bmi.R:** This R script takes a text file containing individual ID numbers and corresponding *BMI* values for each and randomly splits the file into two separate data sets. The seed of R's random number generator is set to `1998`.  

- ***Input:*** bmi4plink.txt 
- ***Output:*** bmi_gwasA.txt, bmi_gwasB.txt

&nbsp;

#### Step 2: Performing two GWASs  


**parallelA_bmi.sh:** This shell script allows us to run the commands contained in *bmiA.sh* for each chromosome in parallel. 

- ***Input:*** bmiA.sh


**parallelB_bmi.sh:** This shell script allows us to run the commands contained in *bmiB.sh* for each chromosome in parallel. 

- ***Input:*** bmiB.sh


**bmiA.sh:** This shell script first uses PLINK 2.0 to produce \*_qcd_bmiA.bim, \*_qcd_bmiA.bed, \*_qcd_bmiA.fam files, with the phenotype data being specified in *bmi_gwasA.txt*. A *linear* model is then fitted for every variant in which 8 principal components are included as well as age.

- ***Input:*** \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, bmi_gwasA.txt, covariates.txt
- ***Output:*** bmiA_res_\*.PHENO1.glm.linear for each chromosome

**bmiB.sh:** This shell script first uses PLINK 2.0 to produce \*_qcd_bmiB.bim, \*_qcd_bmiB.bed, \*_qcd_bmiB.fam files, with the phenotype data being specified in *bmi_gwasB.txt*. A *linear* model is then fitted for every variant in which 8 principal components are included as well as age.

- ***Input:*** \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, bmi_gwasB.txt, covariates.txt
- ***Output:*** bmiB_res_\*.PHENO1.glm.linear for each chromosome

&nbsp;

#### Step 3: Combining results 

**summary_data_bmi.R:** This R script is designed to combine all summary statistics obtained from performing our two *BMI* GWASs together in a suitable format, in which one GWAS is considered to be the discovery GWAS and the other the replication GWAS. The output is a text file which contains information regarding chromosome, position, rsID, effect size and corresponding standard error from the discovery GWAS as well as effect size from the replication GWAS for each variant. 

**summary_data_bmi_1.sh:** This shell script specifies that the files bmiA_res_\*.PHENO1.glm.linear for chromosomes 1 to 22 are to be used as the discovery GWAS results and the files bmiB_res_\*.PHENO1.glm.linear as the replication GWAS results for the R script summary_data_bmi.R.

- ***Input:*** summary_data_bmi.R, bmiA_res_\*.PHENO1.glm.linear for each chromosome, bmiB_res_\*.PHENO1.glm.linear for each chromosome
- ***Output:*** summary_data_bmi_1.txt

**summary_data_bmi_2.sh:** This shell script specifies that the files bmiB_res_\*.PHENO1.glm.linear for chromosomes 1 to 22 are to be used as the discovery GWAS results and the files *bmiA*_res_\*.PHENO1.glm.linear as the replication GWAS results for the R script summary_data_bmi.R.

- ***Input:*** summary_data_bmi.R, bmiB_res_\*.PHENO1.glm.linear for each chromosome, bmiA_res_\*.PHENO1.glm.linear for each chromosome
- ***Output:*** summary_data_bmi_2.txt

&nbsp;

#### Step 4: Obtaining pruned set of SNPs

A pruned set of SNPs was also obtained. Pruning occurred by first calculating LD between each pair of SNPs in a window of 50 SNPs. If an LD value greater than 0.5 was observed, then one SNP out of this pair was removed. The window was shifted 5 SNPs forward and the process was repeated. 1,589,295 SNPs remained after this procedure, a data set about 20% of the size of the original. 

The following code was run: `sbatch linkdis_bmiA.sh` 

**join_pruned.R:** This R script is designed to combine the lists of pruned SNPs from each chromosome together to obtain one single list of pruned SNPs. 

**linkdis_bmiA.sh:** This shell script performs pruning using the command ‘--indep-pairwise 50 5 0.5’ to obtain a full list of pruned SNPs.

- ***Input:*** \*_qcd_bmiA.bim, \*_qcd_bmiA.bed, \*_qcd_bmiA.fam files for each chromosome
- ***Output:*** pruned_SNPs_bmi_1.txt

&nbsp;

### Application of Winner's Curse methods 

#### Evaluating methods

**winnerscurse_realdata.R:** This script uses 6 sets of GWAS summary statistics, two related to each trait, which were generated in the previous section. This script's primary purpose is to evaluate a number of winner's curse correction methods by comparing their performance at two significance thresholds, 5 <span>&#215;</span> 10<sup>-8</sup> and 5 <span>&#215;</span> 10<sup>-4</sup>, for each data set, using both the estimated MSE among significant SNPs and the average bias over all significant SNPs. It also performs an initial exploration of the data sets by computing the number of significant SNPs at each thresholds as well as proportions that indicate the extent of winner's curse.

- ***Input:*** summary_data_bmi_1.txt, summary_data_bmi_2.txt, summary_data_T2D_1.txt, summary_data_T2D_2.txt, summary_data_height_1.txt, summary_data_height_2.txt
- ***Output:*** realdata_Table1.txt, mse_5e_8.txt, mse_5e_4.txt, bias_5e_8_positive.txt, bias_5e_8_negative.txt, bias_5e_4_positive.txt, bias_5e_4_negative.txt

&nbsp;

#### Evaluating methods using pruned data sets

**realdata_pruned.R:** This script is very similar to the above script, 'winnerscurse_realdata.R'. However, in this script, pruned versions of the six real data sets are first obtained. Following this, the winner's curse correction methods are then evaluated at two significance thresholds, 5 <span>&#215;</span> 10<sup>-8</sup> and 5 <span>&#215;</span> 10<sup>-4</sup>, for each pruned data set, using the estimated MSE among significant SNPs.

- ***Input:*** pruned_SNPs_bmi_1.txt, summary_data_bmi_1.txt, summary_data_bmi_2.txt, summary_data_T2D_1.txt, summary_data_T2D_2.txt, summary_data_height_1.txt, summary_data_height_2.txt
- ***Output:*** pruned_mse_5e_8.txt, pruned_mse_5e_4.txt

&nbsp;

#### Investigating number of independent signals required

**ind_signals.R:** This script accompanies our attempt at determining the number of independent signals required to ensure appropriate performance of the winner's curse correction methods. Firstly, the estimated MSE is computed for the BMI data sets at thresholds 5 <span>&#215;</span> 10<sup>-10</sup>, 5 <span>&#215;</span> 10<sup>-12</sup> and 5 <span>&#215;</span> 10<sup>-14</sup>, and for the height data sets at thresholds 5 <span>&#215;</span> 10<sup>-32</sup>, 5 <span>&#215;</span> 10<sup>-34</sup> and 5 <span>&#215;</span> 10<sup>-36</sup>. These results are illustrated. Manhattan plots for each data set are also produced. 

- ***Input:*** summary_data_bmi_1.txt, summary_data_bmi_2.txt, summary_data_height_1.txt, summary_data_height_2.txt
- ***Output:*** S23_Fig.tiff, S24_Fig.tiff, S25_Fig.tiff, S26_Fig.tiff, S27_Fig.tiff, S28_Fig.tiff 

&nbsp;

## Figures

The two scripts below were used to produce all figures in the manuscript, as well as supplementary figures, and can be found in the `figures` folder.

**sims_figs.R:** This script was used to generate the following figures: Fig1.tiff, S1_Fig.tiff - S15_Fig.tiff, S29_Fig.tiff

**realdata_figs.R:** This script was used to generate the following figures: Fig2.tiff, Fig3.tiff, S16_Fig.tiff - S22_Fig.tiff
