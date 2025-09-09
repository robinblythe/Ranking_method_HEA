# Supplementary file - sensitivity analysis
options(scipen = 100, digits = 3)
library(tidyverse)
library(patchwork)

g_colours <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")

results <- bind_rows(readRDS(file = "sim_results.RDS")) |>
  mutate(Mean_rank = ifelse(is.nan(Mean_rank), Inf, Mean_rank),
         Mean_risk = ifelse(is.nan(Mean_risk), 0, Mean_risk)) |>
  group_by(Strategy, auc_model, Prevalence, Miscalibration) |>
  summarise(Mean_rank = median(Mean_rank),
            Mean_rank_low = quantile(Mean_rank, 0.25),
            Mean_rank_high = quantile(Mean_rank, 0.75),
            Mean_risk = median(Mean_risk),
            Mean_risk_low = quantile(Mean_risk, 0.25),
            Mean_risk_high = quantile(Mean_risk, 0.75),
            PPV_median = median(PPV),
            PPV_low = quantile(PPV, 0.25),
            PPV_high = quantile(PPV, 0.75),
            Sens_median = median(sensitivity),
            Sens_low = quantile(sensitivity, 0.25),
            Sens_high = quantile(sensitivity, 0.75))

p <- results |>
  ggplot(aes(x = auc_model))

(p +
    geom_line(aes(y = Mean_rank, colour = Strategy, linetype = Miscalibration), linewidth = 1.2) +
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
    labs(y = "Mean rank of true positives")) /
  (p +
     geom_line(aes(y = Mean_risk, colour = Strategy, linetype = Miscalibration), linewidth = 1.2) +
     facet_wrap(vars(Prevalence), labeller = label_both) +
     theme_bw() +
     theme(panel.grid.minor = element_blank(), 
           axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           legend.position = "none") +
     scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
     scale_colour_manual(values = g_colours) +
     scale_fill_manual(values = g_colours) +
     labs(y = "Mean predicted risk")) /
  (p +
     geom_line(aes(y = PPV_median, colour = Strategy, linetype = Miscalibration), linewidth = 1.2) +
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
     geom_line(aes(y = Sens_median, colour = Strategy, linetype = Miscalibration), linewidth = 1.2) +
     facet_wrap(vars(Prevalence), labeller = label_both) +
     theme_bw() +
     labs(x = "Model AUC",
          y = "Sensitivity") +
     scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
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
    scale_x_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, 0.1)) +
    scale_y_continuous(labels = scales::dollar_format(big.mark = ",")) +
    scale_colour_manual(values = g_colours) +
    scale_fill_manual(values = g_colours) +
    labs(y = "False positive cost (SGD)") +
    ggtitle("Risks overestimated")) /
  (p2 +
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
  (p2 +
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

ggsave(filename = "Figure S2.jpg", width = 8, height = 8)

