obtain_class <- function(data) {

  do.call(rbind, data) |>
    as_tibble() |>
    mutate(predicted = predict(fit, newdata = pick(x), type = "response"),
           class_val_opt = ifelse(predicted >= cutpoint_nmb, 1, 0),
           class_nne = ifelse(predicted >= cutpoint_nne, 1, 0),
           class_youden = ifelse(predicted >= cutpoint_youden, 1, 0)
           )
}

obtain_sample <- function(data, threshold_class, sample){

  subset(data, get(threshold_class) == 1) |>
    group_by(auc, p0) |>
    slice_sample(n = sample) |>
    ungroup() |>
    mutate(Method = gsub("class_", "", threshold_class)) |>
    select(Method, auc, p0, predicted, actual)
  
}

simulator <- function(data, n_samples, n_eval, max_iter){
  results_out <- list()
  for (i in 1:max_iter) {
    
    # Sample from data and take all patients with a score above t1 but below t3
    df = df_ed[sample(nrow(df_ed), n_samples),] |>
      filter(t1 == 1 & t3 == 0) |>
      arrange(desc(pred_risk))
    
    FP = rgamma(nrow(df), shape = 110.314, scale = 0.172) * (3.19 * 0.85)
    df$Cost = FP
    
    df_rand = df |>
      slice_sample(n = n_eval) |>
      mutate(
        iter = i,
        Method = "Unranked")
    
    df_rank = df |>
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
  