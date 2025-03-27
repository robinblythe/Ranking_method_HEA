### Functions for simulation study
# Helper function to get classification
obtain_class <- function(data, fit) {
  do.call(rbind, data) |>
    as_tibble() |>
    mutate(
      predicted = predict(fit, newdata = pick(x), type = "response"),
      class_val_opt = ifelse(predicted >= cutpoint_nmb, 1, 0),
      class_nne = ifelse(predicted >= cutpoint_nne, 1, 0),
      class_youden = ifelse(predicted >= cutpoint_youden, 1, 0)
    )
}

# Helper function to get everyone in a class and sample from that group randomly
obtain_sample <- function(data, threshold_class, sample) {
  subset(data, get(threshold_class) == 1) |>
    group_by(auc, p0) |>
    slice_sample(n = sample) |>
    ungroup() |>
    mutate(Method = gsub("class_", "", threshold_class)) |>
    select(Method, auc, p0, predicted, actual)
}

# Helper function for ROC-based threshold
# Obtain a cutpoint based on number needed to evaluate (NNE)
# Equal to maximum tolerable number of false positives per true positive
get_nne_threshold <- function(predictor, response, nne){
  roc_curve <- roc(predictor = predictor, response = response)
  ppv <- coords(roc_curve, x = "all", input = "threshold", ret = c("threshold", "ppv"))
  ppv$nne <- 1 / ppv$ppv
  
  return(min(ppv$threshold[ppv$nne == nne], na.rm = T))
}

# Main simulation function to repeat the analysis across different event rates, aucs, sample sizes
run_sims <- function(event_rate, auc, samp_size_multi, niter, n_test, n_eval) {
  # Get sample size requirements (once per input combination)
  sampsize = pmsampsize::pmsampsize(type = "b", parameters = 1, prevalence = event_rate, cstatistic = auc)
  
  results <- list()
  # Run for loop 
  for (i in 1:niter) {
    # Simulate model training data
    train = predictNMB::get_sample(
      auc = auc,
      n_samples = sampsize$sample_size * samp_size_multi,
      prevalence = event_rate,
      min_events = 0
    )
    # Fit model
    fit = glm(actual ~ x, data = train, family = binomial())
    
    # Obtain predictions from fitted model to select thresholds
    train$predicted = predict(fit, type = "response")
    
    # Run predictNMB simulation
    nmb_sim = do_nmb_sim(
      sample_size = nrow(train),
      n_sims = 50,
      n_valid = 10000,
      sim_auc = auc,
      event_rate = sum(train$actual == 1)/nrow(train),
      fx_nmb_training = fx_nmb,
      fx_nmb_evaluation = fx_nmb
    )
    
    # Extract youden index cutpoints
    cutpoint_youden = median(nmb_sim$df_thresholds$youden)
    cutpoint_nmb = median(nmb_sim$df_thresholds$value_optimising)
    
    # Extract Number Needed to Evaluate (NNE) cutpoint from ROC curve
    cutpoint_nne = get_nne_threshold(predictor = train$predicted, response = train$actual, nne = nne)
    
    # Simulate model validation data - e.g., a 1000 bed hospital monitoring deteriorating patients
    test = get_sample(auc, n_test, event_rate)
    test$predicted = predict(fit, type = "response", newdata = test)
    FP = rgamma(1, shape = 110.314, scale = 0.172) * 3.19 +
      (event_rate * (1 - rnorm(1, 0.910, 0.036)) * rnorm(1, 14134, 686))
    sample_youden = subset(test, predicted >= cutpoint_youden) |> slice_sample(n = n_eval)
    sample_nmb = subset(test, predicted >= cutpoint_nmb) |> slice_sample(n = n_eval)
    sample_nne = subset(test, predicted >= cutpoint_nne) |> slice_sample(n = n_eval)
    sample_rank = test |> arrange(desc(predicted)) |> slice_head(n = n_eval)
    
    results[[i]] <- tibble(
      iter = i,
      strategy = c("Youden", "NMB", "NNE", "Rank"),
      FP_cost = c(
        sum(sample_youden$actual == 0) * FP,
        sum(sample_nmb$actual == 0) * FP,
        sum(sample_nne$actual == 0) * FP,
        sum(sample_rank$actual == 0) * FP
      ),
      auc = auc,
      event_rate = event_rate,
      pr_required_sampsize = samp_size_multi
    )
  }
  
  return(bind_rows(results))
}




### Functions for real data example (SERP)
# Function to sample from data and return ranked/unranked datasets
simulator <- function(data, n_samples, n_eval, max_iter) {
  results_out <- list()
  for (i in 1:max_iter) {
    # Sample from data and take all patients with a score above t1 but below t3
    df <- df_ed[sample(nrow(df_ed), n_samples), ] |>
      filter(t1 == 1 & t3 == 0) |>
      arrange(desc(pred_risk))

    FP <- rgamma(nrow(df), shape = 110.314, scale = 0.172) * (3.19 * 0.85)
    df$Cost <- FP

    df_rand <- df |>
      slice_sample(n = n_eval) |>
      mutate(
        iter = i,
        Method = "Unranked"
      )

    df_rank <- df |>
      slice_head(n = n_eval) |>
      mutate(
        iter = i,
        Method = "Ranked"
      )

    results_out[[i]] <- do.call(rbind, list(df_rand, df_rank))
  }

  serp_results <- do.call(rbind, results_out)
  serp_results$TP <- ifelse(serp_results$class_t1 == "TP", 1, 0)
  serp_results$Cost <- if_else(serp_results$TP == 1, 0, serp_results$Cost)

  return(serp_results)
}
