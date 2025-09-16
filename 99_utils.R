# Helper function for ROC-based threshold
# Obtain a cutpoint based on number needed to evaluate (NNE)
# Equal to maximum tolerable number of false positives per true positive
get_nne_threshold <- function(predictor, response, nne) {
  roc_curve <- pROC::roc(predictor = predictor, response = response)
  ppv <- pROC::coords(
    roc_curve,
    x = "all", input = "threshold", ret = c("threshold", "ppv")
  )
  ppv$nne <- 1 / ppv$ppv
  cutpoint <- ppv$threshold[which.min(abs(ppv$nne - nne))]
  
  if_else(is.infinite(cutpoint), 0, cutpoint)
}

# Main simulation function to repeat the analysis across different event rates, aucs, sample sizes
run_sims <- function(event_rate, auc, miscalibration, niter, n_test, n_eval, seed) {
  
  # Set up health economic helper function for predictNMB value-optimising cutpoint
  wtp <- 45000
  nmb = predictNMB::get_nmb_sampler(
    # Cost of ICU admission
    outcome_cost = 14345,
    # Willingness to pay per QALY
    wtp = wtp,
    # QALYs lost due to deterioration event
    qalys_lost = 0.03,
    # Cost of an evaluation = (Clinician time cost * duration of evaluation) + (Opportunity cost = chance of successful intervention * outcome cost * underlying p0)
    high_risk_group_treatment_cost = (1.50 * 19) + ((1 - 0.910) * 14345 * event_rate),
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
      n_samples = sampsize$sample_size,
      prevalence = event_rate,
      min_events = ceiling(sampsize$events)
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
    test$predicted_calibrated = predict(fit, type = "response", newdata = test)
    
    # Induce miscalibration if specified
    df_position <- test |> 
      mutate(predicted = 
               case_when(
                 miscalibration == "underestimates" ~ predicted_calibrated^alpha,
                 miscalibration == "overestimates" ~ 1 - (1 - predicted_calibrated)^alpha,
                 .default = predicted_calibrated
                 )
             ) |>
      # Arrange by rank
      arrange(desc(predicted)) |>
      mutate(
        Youden_pos = predicted >= cutpoint_youden,
        NMB_pos = predicted >= cutpoint_nmb,
        NNE_pos = predicted >= cutpoint_nne,
        Rank_pos = row_number() <= n_eval
      ) |>
      pivot_longer(
        ends_with("pos"),
        names_to = "Method",
        values_to = "above_threshold",
        names_transform = ~ str_remove(.x, "_pos")
      )
  
    results[[i]] = df_position |>
      filter(above_threshold) |>
      group_by(Method) |>
      mutate(
        rank_rand = sample(1:n(), n(), replace = F)
      ) |>
      filter(
        case_when(Method != "Rank" ~ rank_rand <= n_eval,
        .default = above_threshold
        )) |>
      summarise(
        ppv = sum(actual) / n(),
        sensitivity = sum(actual) / sum(test$actual),
        mean_rank = if (unique(Method) == "Rank") {
          mean(row_number()[actual == 1])
        } else {
          mean(rank_rand[actual == 1])
        },
        mean_risk = mean(predicted_calibrated),
        .groups = "drop"
      ) |>
      mutate(
        auc_model = auc,
        Prevalence = event_rate,
        Miscalibration = miscalibration,
        iter = i
      )
  }
  return(bind_rows(results))
}



# Applied simulation function (SERP)
run_serp <- function(niter, cutpoint, seed){
  results <- list()
  set.seed(seed)
  for (i in 1:niter){
    
    # Study costs originally reported in AUD at different time points
    # Converted from base year to 2025 (Cost * (1.03)^(2025 - year))
    # Converted from AUD to SGD (1 AUD = 0.85 SGD as of 9 September 2025)
    deterioration_cost = rnorm(1, 14345, 696) # see Curtis et al
    
    FP = rgamma(1, shape = 110.314, scale = 0.172) * 1.50 + # cost of clinical time - see Mudiyaselange et al who report 2016 costs
      (event_rate * (1 - rnorm(1, 0.910, 0.036)) * deterioration_cost) # opportunity cost
    FN = deterioration_cost
    df_sample = df_ed[sample(nrow(df_ed), 150),]
    df_positive <- df_sample |>
      filter(pred_risk >= cutpoint) |>
      arrange(desc(pred_risk)) |>
      mutate(rank = row_number(),
             rank_random = sample(1:n(), n(), replace = FALSE),
             top_25_ranked = ifelse(rank <= 0.25 * n(), TRUE, FALSE),
             top_25_unranked = ifelse(rank_random <= 0.25 * n(), TRUE, FALSE),
             top_50_ranked = ifelse(rank <= 0.5 * n(), TRUE, FALSE),
             top_50_unranked = ifelse(rank_random <= 0.5 * n(), TRUE, FALSE),
             top_75_ranked = ifelse(rank <= 0.75 * n(), TRUE, FALSE),
             top_75_unranked = ifelse(rank_random <= 0.75 * n(), TRUE, FALSE)) |>
      pivot_longer(
        starts_with("top"),
        names_to = "constraint",
        values_to = "above_threshold",
        names_transform = ~ str_remove(.x, "_limit")
      ) |>
      filter(above_threshold)
    
    results[[i]] <- df_positive |>
      summarise(
        .by = constraint,
        PPV = sum(outcome_died_30d)/n(),
        Sensitivity = sum(outcome_died_30d)/sum(df_sample$outcome_died_30d),
        FP_cost = sum(outcome_died_30d == 0) * FP[1],
        Misclassification_cost = sum(outcome_died_30d == 0) * FP[1] +
          (sum(df_sample$outcome_died_30d == 1) - sum(outcome_died_30d == 1)) * FN[1]
      )
  }
  
  return(bind_rows(results))
  
}
