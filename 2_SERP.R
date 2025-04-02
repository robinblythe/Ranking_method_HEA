options(scipen = 100, digits = 5)

library(ggplot2)
library(patchwork)

# Load ED data
df_ed <- readRDS("C:/Users/blythe/NUS Dropbox/EDData/ED2A/ED2A_Blythe Robin Daniel/dat_scored.RDS")

# Assess external validity for 30d mortality (purpose of model) and ICU admission (not purpose of model)
preds <- subset(df_ed, select = c("outcome_died_30d", "outcome_icu", "pred_risk"))
preds$pr_30d <- predict(glm(outcome_died_30d ~ pred_risk, family = binomial(), data = preds), type = "response")
preds$pr_icu <- predict(glm(outcome_icu ~ pred_risk, family = binomial(), data = preds), type = "response")


# 30d mortality
p <- preds |> ggplot()
(p +
  geom_histogram(aes(x = pr_30d, y = after_stat(count)/sum(after_stat(count))), alpha = 0.2, colour = "#999999") +
  geom_smooth(aes(x = pr_30d, y = outcome_died_30d)) +
  geom_abline() +
  theme_bw() +
  ylab("Observed frequency") +
  xlab("Predicted probability") +
  geom_text(aes(x = 0.3, y = 0.9), label = paste0("Model AUC (30d mortality) = ", 
                                                   round(auc(predictor = preds$pr_30d, response = preds$outcome_died_30d)[[1]], 2))) +
  geom_text(aes(x = 0.4, y = 0.65), label = "Model underestimates risks") +
  geom_text(aes(x = 0.6, y = 0.25), label = "Model overestimates risks")) /
(p +
   geom_histogram(aes(x = pr_icu, y = after_stat(count)/sum(after_stat(count))), alpha = 0.2, colour = "#999999") +
   geom_smooth(aes(x = pr_icu, y = outcome_icu), colour = "red") +
   geom_abline() +
   theme_bw() +
   ylab("Observed frequency") +
   xlab("Predicted probability") +
   geom_text(aes(x = 0.3, y = 0.9), label = paste0("Model AUC (ICU admission) = ", 
                                                   round(auc(predictor = preds$pr_icu, response = preds$outcome_icu)[[1]], 2))) +
   geom_text(aes(x = 0.4, y = 0.65), label = "Model underestimates risks") +
   geom_text(aes(x = 0.6, y = 0.25), label = "Model overestimates risks")) +
  plot_annotation(tag_levels = 'A')
     
ggsave(filename = "Figure 2.jpg", height = 8, width = 6)


# Threshold value of 27 as per https://www.medrxiv.org/content/10.1101/2021.02.09.21251397v1.full
cutpoint <- 27