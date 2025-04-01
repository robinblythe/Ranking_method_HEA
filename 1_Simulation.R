# Set up study
options(scipen = 999, digits = 3)
library(predictNMB)
library(parallel)
library(pROC)
library(tidyverse)

source("./99_utils.R")

# Select number needed to evaluate for the NNE threshold
nne <- 14

# Set up parallel computing
# cl = makeCluster(detectCores() - 2)

# Simulation values:
n_test <- 1000
n_eval <- 40
niter <- 10000

test <- run_sims(event_rate = 0.05, auc = 0.85, samp_size_multi = 0.8, niter = niter, n_test = n_test, n_eval = n_eval)


