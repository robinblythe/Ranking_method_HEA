# Supplementary file - sensitivity analysis
options(scipen = 100, digits = 3)
library(tidyverse)
library(patchwork)

g_colours <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")

results <- bind_rows(readRDS(file = "sim_results.RDS")) |>
  group_by(Strategy, auc_model, Prevalence, pr_required_sampsize) |>
  summarise(FP_cost_median = median(FP_cost),
            FP_cost_low = quantile(FP_cost, 0.25),
            FP_cost_high = quantile(FP_cost, 0.75),
            PPV_median = median(PPV),
            PPV_low = quantile(PPV, 0.25),
            PPV_high = quantile(PPV, 0.75),
            Sens_median = median(sensitivity),
            Sens_low = quantile(sensitivity, 0.25),
            Sens_high = quantile(sensitivity, 0.75))

p1 <- results |>
  filter(pr_required_sampsize == 0.8) |>
  ggplot(aes(x = auc_model))

p2 <- results |>
  filter(pr_required_sampsize == 1.2) |>
  ggplot(aes(x = auc_model))

(p1 +
    geom_line(aes(y = FP_cost_median, colour = Strategy), linewidth = 1.2) +
    geom_ribbon(aes(ymin = FP_cost_low, ymax = FP_cost_high, fill = Strategy), alpha = 0.2) +
    facet_wrap(vars(Prevalence), labeller = label_both) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none") +
    scale_y_continuous(labels = scales::dollar_format(big.mark = ",")) +
    scale_colour_manual(values = g_colours) +
    scale_fill_manual(values = g_colours) +
    labs(y = "False positive cost (SGD)") +
    ggtitle("Underpowered: Sample size 80% of required sample")) /
  (p1 +
     geom_line(aes(y = PPV_median, colour = Strategy), linewidth = 1.2) +
     geom_ribbon(aes(ymin = PPV_low, ymax = PPV_high, fill = Strategy), alpha = 0.2) +
     facet_wrap(vars(Prevalence), labeller = label_both) +
     theme_bw() +
     theme(panel.grid.minor = element_blank(), 
           axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = "none") +
     scale_colour_manual(values = g_colours) +
     scale_fill_manual(values = g_colours) +
     labs(y = "Positive Predictive Value")) /
  (p1 +
     geom_line(aes(y = Sens_median, colour = Strategy), linewidth = 1.2) +
     geom_ribbon(aes(ymin = Sens_low, ymax = Sens_high, fill = Strategy), alpha = 0.2) +
     facet_wrap(vars(Prevalence), labeller = label_both) +
     theme_bw() +
     labs(x = "Model AUC",
          y = "Sensitivity") +
     scale_colour_manual(values = g_colours) +
     scale_fill_manual(values = g_colours) +
     theme(panel.grid.minor = element_blank(),
           legend.position = "bottom"))

ggsave(filename = "Figure S1.jpg", width = 8, height = 8)


(p2 +
    geom_line(aes(y = FP_cost_median, colour = Strategy), linewidth = 1.2) +
    geom_ribbon(aes(ymin = FP_cost_low, ymax = FP_cost_high, fill = Strategy), alpha = 0.2) +
    facet_wrap(vars(Prevalence), labeller = label_both) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none") +
    scale_y_continuous(labels = scales::dollar_format(big.mark = ",")) +
    scale_colour_manual(values = g_colours) +
    scale_fill_manual(values = g_colours) +
    labs(y = "False positive cost (SGD)") +
    ggtitle("Overpowered: Sample size 120% of required sample")) /
  (p2 +
     geom_line(aes(y = PPV_median, colour = Strategy), linewidth = 1.2) +
     geom_ribbon(aes(ymin = PPV_low, ymax = PPV_high, fill = Strategy), alpha = 0.2) +
     facet_wrap(vars(Prevalence), labeller = label_both) +
     theme_bw() +
     theme(panel.grid.minor = element_blank(), 
           axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = "none") +
     scale_colour_manual(values = g_colours) +
     scale_fill_manual(values = g_colours) +
     labs(y = "Positive Predictive Value")) /
  (p2 +
     geom_line(aes(y = Sens_median, colour = Strategy), linewidth = 1.2) +
     geom_ribbon(aes(ymin = Sens_low, ymax = Sens_high, fill = Strategy), alpha = 0.2) +
     facet_wrap(vars(Prevalence), labeller = label_both) +
     theme_bw() +
     labs(x = "Model AUC",
          y = "Sensitivity") +
     scale_colour_manual(values = g_colours) +
     scale_fill_manual(values = g_colours) +
     theme(panel.grid.minor = element_blank(),
           legend.position = "bottom"))

ggsave(filename = "Figure S2.jpg", width = 8, height = 8)

