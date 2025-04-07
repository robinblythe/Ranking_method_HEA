### Functions for simulation study
library(tidyverse)

# Helper function for ROC-based threshold
# Obtain a cutpoint based on number needed to evaluate (NNE)
# Equal to maximum tolerable number of false positives per true positive
get_nne_threshold <- function(predictor, response, nne){
  roc_curve = pROC::roc(predictor = predictor, response = response)
  ppv = pROC::coords(roc_curve, x = "all", input = "threshold", ret = c("threshold", "ppv"))
  ppv$nne = 1 / ppv$ppv
  cutpoint = min(ppv$threshold[ppv$nne == nne], na.rm = T) |> suppressWarnings()
  
  return(ifelse(is.infinite(cutpoint), 0, cutpoint))
}

# Main simulation function to repeat the analysis across different event rates, aucs, sample sizes
run_sims <- function(event_rate, auc, samp_size_multi, niter, n_test, n_eval, seed) {
  
  # Set up health economic helper function for predictNMB value-optimising cutpoint
  wtp <- 45000
  nmb = predictNMB::get_nmb_sampler(
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
  
  # Get sample size requirements (once per input combination)
  sampsize = pmsampsize::pmsampsize(type = "b", parameters = 2, prevalence = event_rate, cstatistic = auc)
  
  results = list()
  set.seed(seed)
  # Run for loop 
  for (i in 1:niter) {
    
    # Simulate model training data
    train = predictNMB::get_sample(
      auc = auc,
      n_samples = sampsize$sample_size * samp_size_multi,
      prevalence = event_rate,
      min_events = ceiling(sampsize$events) * samp_size_multi
    )
    
    # Fit model
    fit = glm(actual ~ x, data = train, family = binomial())
    
    # Obtain predictions from fitted model to select thresholds
    train$predicted = predict(fit, type = "response")
    
    # Run predictNMB simulation
    thresholds = predictNMB::get_thresholds(predicted = train$predicted, 
                                            actual = train$actual,
                                            nmb = c("TP" = nmb()[["TP"]],
                                                    "TN" = nmb()[["TN"]],
                                                    "FP" = nmb()[["FP"]],
                                                    "FN" = nmb()[["FN"]]))
    
    # Extract youden index cutpoints
    cutpoint_nmb = thresholds[["value_optimising"]]
    cutpoint_youden = thresholds[["youden"]]
    
    # Extract Number Needed to Evaluate (NNE) cutpoint from ROC curve
    cutpoint_nne = get_nne_threshold(predictor = train$predicted, response = train$actual, nne = nne) |> suppressMessages()
    
    # Simulate model validation data - e.g., a 1000 bed hospital monitoring deteriorating patients
    test = predictNMB::get_sample(auc, n_test, event_rate)
    test$predicted = predict(fit, type = "response", newdata = test)
    FP = rgamma(1, shape = 110.314, scale = 0.172) * 3.19 +
      (event_rate * (1 - rnorm(1, 0.910, 0.036)) * rnorm(1, 14134, 686))
    sample_youden = subset(test, predicted >= cutpoint_youden) |> slice_sample(n = n_eval)
    sample_nmb = subset(test, predicted >= cutpoint_nmb) |> slice_sample(n = n_eval)
    sample_nne = subset(test, predicted >= cutpoint_nne) |> slice_sample(n = n_eval)
    sample_rank = test |> arrange(desc(predicted)) |> slice_head(n = n_eval)
    
    results[[i]] <- tibble(
      iter = i,
      Strategy = c("Youden", "NMB", "NNE", "Rank"),
      PPV = c(
        sum(sample_youden$actual)/n_eval,
        sum(sample_nmb$actual)/n_eval,
        sum(sample_nne$actual)/n_eval,
        sum(sample_rank$actual)/n_eval
      ),
      sensitivity = c(
        sum(sample_youden$actual)/sum(test$actual),
        sum(sample_nmb$actual)/sum(test$actual),
        sum(sample_nne$actual)/sum(test$actual),
        sum(sample_rank$actual)/sum(test$actual)
      ),
      FP_cost = c(
        sum(sample_youden$actual == 0) * FP,
        sum(sample_nmb$actual == 0) * FP,
        sum(sample_nne$actual == 0) * FP,
        sum(sample_rank$actual == 0) * FP
      ),
      auc_model = auc,
      Prevalence = event_rate,
      pr_required_sampsize = samp_size_multi
    )
  }
  
  return(bind_rows(results))
}



# Applied simulation function (SERP)
run_serp <- function(niter, cutpoint, seed){
  results <- list()
  set.seed(seed)
  for (i in 1:niter){
    FP = rgamma(1, shape = 110.314, scale = 0.172) * 3.19 +
      (event_rate * (1 - rnorm(1, 0.910, 0.036)) * rnorm(1, 14134, 686))
    df_sample = df_ed[sample(nrow(df_ed), 150),]
    df_positive = subset(df_sample, pred_risk >= cutpoint) |> arrange(desc(pred_risk))
    
    sample_25_rank = df_positive |> slice_head(n = floor(nrow(df_positive) * 0.25))
    sample_50_rank = df_positive |> slice_head(n = floor(nrow(df_positive) * 0.50))
    sample_75_rank = df_positive |> slice_head(n = floor(nrow(df_positive) * 0.75))
    sample_25_threshold = df_positive |> slice_sample(n = floor(nrow(df_positive) * 0.25))
    sample_50_threshold = df_positive |> slice_sample(n = floor(nrow(df_positive) * 0.50))
    sample_75_threshold = df_positive |> slice_sample(n = floor(nrow(df_positive) * 0.75))
    
    results[[i]] <- tibble(
      iter = i,
      Strategy = c(rep("Ranked", 3), rep("Unranked", 3)),
      Pct_evaluated = rep(c(25, 50, 75), 2),
      PPV = c(
        sum(sample_25_rank$outcome_died_30d)/nrow(sample_25_rank),
        sum(sample_50_rank$outcome_died_30d)/nrow(sample_50_rank),
        sum(sample_75_rank$outcome_died_30d)/nrow(sample_75_rank),
        sum(sample_25_threshold$outcome_died_30d)/nrow(sample_25_threshold),
        sum(sample_50_threshold$outcome_died_30d)/nrow(sample_50_threshold),
        sum(sample_75_threshold$outcome_died_30d)/nrow(sample_75_threshold)
      ),
      Sensitivity = c(
        sum(sample_25_rank$outcome_died_30d)/sum(df_sample$outcome_died_30d),
        sum(sample_50_rank$outcome_died_30d)/sum(df_sample$outcome_died_30d),
        sum(sample_75_rank$outcome_died_30d)/sum(df_sample$outcome_died_30d),
        sum(sample_25_threshold$outcome_died_30d)/sum(df_sample$outcome_died_30d),
        sum(sample_50_threshold$outcome_died_30d)/sum(df_sample$outcome_died_30d),
        sum(sample_75_threshold$outcome_died_30d)/sum(df_sample$outcome_died_30d)
      ),
      FP_cost = c(
        sum(sample_25_rank$outcome_died_30d == 0) * FP,
        sum(sample_50_rank$outcome_died_30d == 0) * FP,
        sum(sample_75_rank$outcome_died_30d == 0) * FP,
        sum(sample_25_threshold$outcome_died_30d == 0) * FP,
        sum(sample_50_threshold$outcome_died_30d == 0) * FP,
        sum(sample_75_threshold$outcome_died_30d == 0) * FP
      )
    )
  }
  
  return(bind_rows(results))
  
}
