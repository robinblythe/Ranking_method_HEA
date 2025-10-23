# Set up study
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
alpha <- 3
seed <- 888

# Set up values for simulations
combs <- expand.grid(
  event_rate = c(0.01, 0.05, 0.1),
  auc = c(0.65, 0.75, 0.85, 0.95),
  miscalibration = c("none", "underestimates", "overestimates")
)

# Run in parallel using futures
plan(multisession, workers = availableCores())

results <- future_map(1:nrow(combs), function(i){
  params = combs[i,]
  run_sims(event_rate = params$event_rate, 
           auc = params$auc, 
           miscalibration = params$miscalibration, 
           niter = niter, n_test = n_test, n_eval = n_eval, seed = seed)
}, .options = furrr_options(seed = TRUE)
)

saveRDS(results, file = "sim_results.rds")

# Miscalibration is explored in the sensitivity analysis
df_summary <- bind_rows(results) |>
  filter(Miscalibration == "none") |>
  select(-Miscalibration) |>
  mutate(Mean_rank = ifelse(is.nan(mean_rank), Inf, mean_rank),
         Mean_risk = ifelse(is.nan(mean_risk), 0, mean_risk)) |>
  group_by(Method, auc_model, Prevalence) |>
  summarise(
    Mean_rank = median(Mean_rank),
    Mean_risk = median(Mean_risk),
    PPV_median = median(ppv),
    PPV_low = quantile(ppv, 0.25),
    PPV_high = quantile(ppv, 0.75),
    Sens_median = median(sensitivity),
    Sens_low = quantile(sensitivity, 0.25),
    Sens_high = quantile(sensitivity, 0.75))

p <- df_summary |>
  ggplot(aes(x = auc_model))

g_colours <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")


(p +
    geom_line(aes(y = Mean_rank, colour = Method), linewidth = 1.2) +
    facet_wrap(vars(Prevalence), labeller = label_both) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none") +
    scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
    scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
    scale_colour_manual(values = g_colours) +
    scale_fill_manual(values = g_colours) +
    labs(y = "Mean rank (true positives)")) /
(p +
    geom_line(aes(y = Mean_risk, colour = Method), linewidth = 1.2) +
    facet_wrap(vars(Prevalence), labeller = label_both) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none") +
    scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
    scale_colour_manual(values = g_colours) +
    scale_fill_manual(values = g_colours) +
    labs(y = "Mean predicted risk (true positives)")) /
(p +
  geom_line(aes(y = PPV_median, colour = Method), linewidth = 1.2) +
  geom_ribbon(aes(ymin = PPV_low, ymax = PPV_high, fill = Method), alpha = 0.2) +
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
  geom_line(aes(y = Sens_median, colour = Method), linewidth = 1.2) +
  geom_ribbon(aes(ymin = Sens_low, ymax = Sens_high, fill = Method), alpha = 0.2) +
  facet_wrap(vars(Prevalence), labeller = label_both) +
  theme_bw() +
  labs(x = "Model AUC",
       y = "Sensitivity") +
  scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
  scale_colour_manual(values = g_colours) +
  scale_fill_manual(values = g_colours) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom"))

ggsave(filename = "./Figures/Figure 2.jpg", height = 10, width = 8)

