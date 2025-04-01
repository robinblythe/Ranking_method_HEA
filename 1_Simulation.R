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
cl = makeCluster(detectCores() - 2)

# Simulation values:
n_test <- 1000
n_eval <- 40

test <- run_sims(event_rate = 0.05, auc = 0.85, samp_size_multi = 0.8, niter = 10, n_test = n_test, n_eval = n_eval)



######################################################################
max_iter <- 100
set.seed(888)
sim_results <- list()

df_sims <- obtain_class(sims, fit)

for (i in 1:max_iter) {
  df_youden <- obtain_sample(df_sims, "class_youden", n_samples)
  df_nne <- obtain_sample(df_sims, "class_nne", n_samples)
  df_val_opt <- obtain_sample(df_sims, "class_val_opt", n_samples)

  df_rank <- df_sims |>
    group_by(auc, p0) |>
    arrange(desc(predicted)) |>
    slice_head(n = n_samples) |>
    mutate(Method = "ranking") |>
    select(Method, auc, p0, predicted, actual)

  iteration <- do.call(rbind, list(df_youden, df_nne, df_val_opt, df_rank)) |>
    group_by(Method, auc, p0) |>
    rowwise() |>
    mutate(
      Outcome = ifelse(actual == 1, "TP", "FP"),
      Cost = ifelse(
        Outcome == "FP",
        # Opportunity cost of a false positive (see methods section)
        rgamma(1, shape = 110.314, scale = 0.172) * 3.19 +
          (p0 * (1 - rnorm(1, 0.910, 0.036)) * rnorm(1, 14134, 686)),
        # No cost for true positives (for now)
        0
      ),
      iter = i
    )

  sim_results[[i]] <- iteration
  remove(df_youden, df_nne, df_val_opt, df_rank, iteration)
}

results <- do.call(rbind, sim_results)
remove(sims)
#########################################################################

