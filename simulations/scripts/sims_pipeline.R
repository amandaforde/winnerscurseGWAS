## WINNER'S CURSE SIMULATION STUDY SCRIPT 1 - PIPELINE:

## This script allows us to run the entire simulation study which evaluates
## and compares various winner's curse correction methods, from beginning to
## end. Note that seeds have been included in each script separately. 

## Load all required packages for simulations.
library(devtools)
devtools::install_github("amandaforde/winnerscurse")
library(winnerscurse)
library(tidyverse)
library(parallel)
library(scam)
library(mgcv)
library(expm)

## Run 'useful_funs.R' in order to define additional functions required for
## the simulations below.
source("simulations/scripts/useful_funs.R")

## Run each set of simulations, taking note of length of time taken for each
## script.
## NOTE: Simulations currently being run on Windows, hence mc.cores=1 in
## mclapply() in below scripts.

start_time <- Sys.time()
source("simulations/scripts/nsig_prop_bias_LD.R")
end_time <- Sys.time()
print("nsig_prop_bias_LD.R complete!")
end_time - start_time

start_time <- Sys.time()
source("simulations/scripts/sims_LD.R")
end_time <- Sys.time()
print("sims_LD.R complete!")
end_time - start_time

start_time <- Sys.time()
source("simulations/scripts/sims_ind_1.R")
end_time <- Sys.time()
print("sims_ind_1.R complete!")
end_time - start_time

start_time <- Sys.time()
source("simulations/scripts/sims_ind_2.R")
end_time <- Sys.time()
print("sims_ind_2.R complete!")
end_time - start_time
