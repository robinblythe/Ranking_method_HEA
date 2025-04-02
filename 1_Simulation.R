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
niter <- 1000

# Set up parallel computing
cl <- makeCluster(detectCores() - 2)

# Set up values for simulations
combs <- expand.grid(
  event_rate = c(0.01, 0.05, 0.1),
  auc = c(0.65, 0.75, 0.85, 0.95),
  samp_size_multi = c(0.8, 1.0, 1.2)
)

results <- mclapply(1:nrow(combs), function(i){
  params = combs[i,]
  run_sims(event_rate = params$event_rate, 
           auc = params$auc, 
           samp_size_multi = params$samp_size_multi, 
           niter = niter, n_test = n_test, n_eval = n_eval, seed = 888)
}
)

saveRDS(results, filename = "sim_results.RDS")

# Do the samp_size_multi as a sensitivity analysis, just stick to the 1x multiplier for main analysis
p <- bind_rows(results) |> 
  filter(pr_required_sampsize == 1) |>
  select(-pr_required_sampsize) |>
  group_by(strategy, auc_model, event_rate_model) |>
  summarise(FP_cost_median = median(FP_cost),
            FP_cost_low = quantile(FP_cost, 0.25),
            FP_cost_high = quantile(FP_cost, 0.75),
            PPV_median = median(PPV),
            PPV_low = quantile(PPV, 0.25),
            PPV_high = quantile(PPV, 0.75),
            Sens_median = median(sensitivity),
            Sens_low = quantile(sensitivity, 0.25),
            Sens_high = quantile(sensitivity, 0.75)) |>
  ggplot(aes(x = auc_model))

p +
  geom_line(aes(y = FP_cost_median, colour = strategy), linewidth = 1.2) +
  geom_ribbon(aes(ymin = FP_cost_low, ymax = FP_cost_high, fill = strategy), alpha = 0.2) +
  facet_wrap(vars(event_rate_model)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = "Figure 1.jpg", height = 8, width = 6)