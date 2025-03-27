# Set up study
options(scipen = 999, digits = 3)
library(predictNMB)
library(parallel)
library(pROC)
library(tidyverse)

source("./99_utils.R")

# Set up health economic helper function for predictNMB value-optimising cutpoint
wtp <- 45000
fx_nmb <- get_nmb_sampler(
  # Cost of ICU admission
  outcome_cost = (14134 * 0.85) * (1.03)^4,
  # Willingness to pay per QALY
  wtp = wtp,
  # QALYs lost due to deterioration event
  qalys_lost = 0.03,
  # Cost of an evaluation = (Clinician time cost * duration of MET) + (Opportunity cost = chance of successful intervention * outcome cost * underlying p0)
  high_risk_group_treatment_cost = (3.19 * 0.85 * 1.03 * 19) + ((1 - 0.910) * (14134 * 0.85) * (1.03)^4 * event_rate),
  # Chance of successful intervention
  high_risk_group_treatment_effect = 1 - 0.910
)

# Select number needed to evaluate
nne <- 14


# Set up parallel computing
cl = makeCluster(detectCores() - 2)




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




# # Notes: run pmsampsize to get the sample size, then multiply it by the samp_size_mult value
#         Also assume we can only ever see 80% of the available patients (can change later)

# for each scenario (combination of event rate and AUC and sample_size_multiplier)

# for each sim (1:10,000):
# sample training data
# fit model
# select threshold (using training data) - only applies to cutpoints, not ranking
# sample hospital scenario data (validation)
# generate predictions for scenario data
# sample false positive cost
# for each strategy (cutpoints and ranking)
# select the n patients recieving treatment
# apply cost to n false positives
# return total cost
# return tibble with one row and one column for each strategy with its cost
