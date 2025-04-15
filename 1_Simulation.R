# Set up study
options(scipen = 100, digits = 3)
library(tidyverse)
library(furrr)
library(patchwork)

source("./99_utils.R")

# Select number needed to evaluate for the NNE threshold
nne <- 14

# Simulation values:
n_test <- 1000
n_eval <- 40
niter <- 10000

# Set up values for simulations
combs <- expand.grid(
  event_rate = c(0.01, 0.05, 0.1),
  auc = c(0.65, 0.75, 0.85, 0.95),
  samp_size_multi = c(0.8, 1.0, 1.2)
)

# Run in parallel using futures
plan(multisession, workers = availableCores())

results <- future_map(1:nrow(combs), function(i){
  params = combs[i,]
  run_sims(event_rate = params$event_rate, 
           auc = params$auc, 
           samp_size_multi = params$samp_size_multi, 
           niter = niter, n_test = n_test, n_eval = n_eval, seed = 888)
}
)

saveRDS(results, file = "sim_results.RDS")

# Do the samp_size_multi as a sensitivity analysis, just stick to the 1x multiplier for main analysis
p <- bind_rows(results) |> 
  filter(pr_required_sampsize == 1) |>
  select(-pr_required_sampsize) |>
  group_by(Strategy, auc_model, Prevalence) |>
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

savings <- with(p$data,
                c("Median savings" = FP_cost_median[Strategy == "Rank" & auc_model == 0.85 & Prevalence == 0.05] - 
                    FP_cost_median[Strategy == "NMB" & auc_model == 0.85 & Prevalence == 0.05],
                  "Savings (Low)" = FP_cost_low[Strategy == "Rank" & auc_model == 0.85 & Prevalence == 0.05] - 
                    FP_cost_low[Strategy == "NMB" & auc_model == 0.85 & Prevalence == 0.05],
                  "Savings (High)" = FP_cost_high[Strategy == "Rank" & auc_model == 0.85 & Prevalence == 0.05] - 
                    FP_cost_high[Strategy == "NMB" & auc_model == 0.85 & Prevalence == 0.05]))

g_colours <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")

(p +
  geom_line(aes(y = FP_cost_median, colour = Strategy), linewidth = 1.2) +
  geom_ribbon(aes(ymin = FP_cost_low, ymax = FP_cost_high, fill = Strategy), alpha = 0.2) +
  facet_wrap(vars(Prevalence), labeller = label_both) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = scales::dollar_format(big.mark = ",")) +
  scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
  scale_colour_manual(values = g_colours) +
  scale_fill_manual(values = g_colours) +
  labs(y = "False positive cost (SGD)")) /
(p +
  geom_line(aes(y = PPV_median, colour = Strategy), linewidth = 1.2) +
  geom_ribbon(aes(ymin = PPV_low, ymax = PPV_high, fill = Strategy), alpha = 0.2) +
  facet_wrap(vars(Prevalence), labeller = label_both) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
  scale_colour_manual(values = g_colours) +
  scale_fill_manual(values = g_colours) +
  labs(y = "Positive Predictive Value")) /
(p +
  geom_line(aes(y = Sens_median, colour = Strategy), linewidth = 1.2) +
  geom_ribbon(aes(ymin = Sens_low, ymax = Sens_high, fill = Strategy), alpha = 0.2) +
  facet_wrap(vars(Prevalence), labeller = label_both) +
  theme_bw() +
  labs(x = "Model AUC",
       y = "Sensitivity") +
  scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
  scale_colour_manual(values = g_colours) +
  scale_fill_manual(values = g_colours) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom"))

ggsave(filename = "Figure 2.jpg", height = 8, width = 8)

