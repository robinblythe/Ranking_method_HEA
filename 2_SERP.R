options(scipen = 100, digits = 5)

library(pROC)
library(tidyverse)
library(patchwork)

# Load ED data
df_ed <- readRDS("C:/Users/blythe/NUS Dropbox/EDData/ED2A/ED2A_Blythe Robin Daniel/dat_scored.RDS")

# Assess external validity for 30d mortality (purpose of model) and ICU admission (not purpose of model)
preds <- subset(df_ed, select = c("outcome_died_30d", "outcome_icu", "pred_risk"))
# Convert to probabilities using logistic regression
preds$pr_30d <- predict(glm(outcome_died_30d ~ pred_risk, family = binomial(), data = preds), type = "response")
preds$pr_icu <- predict(glm(outcome_icu ~ pred_risk, family = binomial(), data = preds), type = "response")

# Plots for external validation
p <- preds |> ggplot()
(p +
  geom_histogram(aes(x = pr_30d, y = after_stat(count)/sum(after_stat(count))), alpha = 0.2, colour = "#999999") +
  geom_smooth(aes(x = pr_30d, y = outcome_died_30d)) +
  geom_abline() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_bw() +
  ylab("Observed frequency") +
  xlab("Predicted probability") +
  geom_text(aes(x = 0.3, y = 0.9), label = paste0("Model AUC (30d mortality) = ", 
                                                   round(auc(predictor = preds$pr_30d, response = preds$outcome_died_30d)[[1]], 2))) +
  geom_text(aes(x = 0.4, y = 0.65), label = "Model underestimates risks") +
  geom_text(aes(x = 0.6, y = 0.25), label = "Model overestimates risks")) #delete this if wanting the combined figure/
(p +
   geom_histogram(aes(x = pr_icu, y = after_stat(count)/sum(after_stat(count))), alpha = 0.2, colour = "#999999") +
   geom_smooth(aes(x = pr_icu, y = outcome_icu), colour = "red") +
   geom_abline() +
   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
   theme_bw() +
   ylab("Observed frequency") +
   xlab("Predicted probability") +
   geom_text(aes(x = 0.3, y = 0.9), label = paste0("Model AUC (ICU admission) = ", 
                                                   round(auc(predictor = preds$pr_icu, response = preds$outcome_icu)[[1]], 2))) +
   geom_text(aes(x = 0.4, y = 0.65), label = "Model underestimates risks") +
   geom_text(aes(x = 0.6, y = 0.25), label = "Model overestimates risks")) +
  plot_annotation(tag_levels = 'A')
     
ggsave(filename = "Figure 3.jpg", height = 5, width = 6)


# Repeat simulation approach
# Threshold value of 27 as per https://www.medrxiv.org/content/10.1101/2021.02.09.21251397v1.full
cutpoint <- 27

# Prevalence of outcome ~ 5%
event_rate <- mean(df_ed$outcome_died_30d)

# Set up and run simulation
niter = 10000

source("./99_utils.R")
results <- run_serp(niter = niter, cutpoint = cutpoint, seed = 888)

saveRDS(results, file = "serp_results.rds")

# Visualise results
p <- results |>
  group_by(Strategy, Pct_evaluated) |>
  summarise(FP_cost_median = median(FP_cost),
            FP_cost_low = quantile(FP_cost, 0.25),
            FP_cost_high = quantile(FP_cost, 0.75),
            PPV_median = median(PPV),
            PPV_low = quantile(PPV, 0.25),
            PPV_high = quantile(PPV, 0.75),
            Sens_median = median(Sensitivity, na.rm = T),
            Sens_low = quantile(Sensitivity, 0.25, na.rm = T),
            Sens_high = quantile(Sensitivity, 0.75, na.rm = T)) |>
  ggplot(aes(x = Pct_evaluated, colour = Strategy, fill = Strategy))

(p +
  geom_line(aes(y = FP_cost_median), linewidth = 1.2) +
  geom_ribbon(aes(ymin = FP_cost_low, ymax = FP_cost_high), alpha = 0.3) +
  scale_colour_manual(values = c("#003D7C", "#EF7C00")) +
  scale_fill_manual(values = c("#003D7C", "#EF7C00")) +
  scale_x_continuous(limits = c(25, 75), breaks = c(25, 50, 75)) +
  theme_bw() +
  labs(y = "False positive cost ($SGD)") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())) /
(p +
   geom_line(aes(y = PPV_median), linewidth = 1.2) +
   geom_ribbon(aes(ymin = PPV_low, ymax = PPV_high), alpha = 0.3) +
   scale_colour_manual(values = c("#003D7C", "#EF7C00")) +
   scale_fill_manual(values = c("#003D7C", "#EF7C00")) +
   scale_x_continuous(limits = c(25, 75), breaks = c(25, 50, 75)) +
   theme_bw() +
   labs(y = "Positive Predictive Value") +
   theme(panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.title.x = element_blank(),
         axis.text.x = element_blank())) /
  (p +
     geom_line(aes(y = Sens_median), linewidth = 1.2) +
     geom_ribbon(aes(ymin = Sens_low, ymax = Sens_high), alpha = 0.3) +
     scale_colour_manual(values = c("#003D7C", "#EF7C00")) +
     scale_fill_manual(values = c("#003D7C", "#EF7C00")) +
     scale_x_continuous(limits = c(25, 75), breaks = c(25, 50, 75)) +
     theme_bw() +
     labs(x = "Percent of patients evaluated",
          y = "Sensitivity") +
     theme(panel.grid.minor = element_blank(),
           legend.position = "bottom")) +
  plot_annotation(tag_levels = 'A')

ggsave("Figure 4.jpg", height = 8, width = 4) 


# Tabulate results
# Note that the false positive costs are for 150 patients. This is equal to 1/900 of the ED's annual flow
summary_table <- p$data |>
  select(Strategy, Pct_evaluated, FP_cost_median, FP_cost_low, FP_cost_high) |>
  pivot_wider(names_from = Strategy, values_from = c(FP_cost_median, FP_cost_low, FP_cost_high)) |>
  mutate(`Savings (annual)` = scales::dollar(round((FP_cost_median_Unranked - FP_cost_median_Ranked) * 900)),
         `95% interval` = paste0("[", scales::dollar((FP_cost_low_Unranked - FP_cost_low_Ranked) * 900), ", ",
                                 scales::dollar((FP_cost_high_Unranked - FP_cost_high_Ranked) * 900), "]")) |>
  select(Pct_evaluated, `Savings (annual)`, `95% interval`) |>
  rename(`% of patients evaluated` = Pct_evaluated)
  
write_csv(summary_table, file = "Annual_savings_SERP.csv")         
         
         
         
         
         
         
         