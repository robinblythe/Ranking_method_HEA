# Set up study
options(scipen = 999, digits = 3)
library(predictNMB)
library(parallel)
library(pROC)
library(tidyverse)

source("./99_utils.R")

# Select number needed to evaluate for the NNE threshold
nne <- 14

# Simulation values:
n_test <- 1000
n_eval <- 40
niter <- 10

# Set up parallel computing
cl <- makeCluster(detectCores() - 2)

# Set up values for simulations
combs <- expand.grid(
  event_rate = c(0.01, 0.05, 0.1),
  auc = c(0.65, 0.75, 0.85, 0.95),
  samp_size_multi = c(0.8, 1.0, 1.2)
)

test <- mclapply(1:nrow(combs), function(i){
  params = combs[i,]
  run_sims(event_rate = params$event_rate, 
           auc = params$auc, 
           samp_size_multi = params$samp_size_multi, 
           niter = niter, n_test = n_test, n_eval = n_eval, seed = 888)
}
)

# Do the samp_size_multi as a sensitivity analysis, just stick to the 1x multiplier for main analysis
p <- bind_rows(test) |> 
  filter(samp_size_multi == 1) |>
  select(-pr_required_sampsize) |>
  ggplot()
