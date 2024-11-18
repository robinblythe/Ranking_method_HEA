obtain_class <- function(df, samples) {

  do.call(rbind, sims) |>
    as_tibble() |>
    mutate(predicted = predict(fit, newdata = pick(x), type = "response"),
           class_val_opt = ifelse(predicted >= cutpoint_nmb, 1, 0),
           class_nne = ifelse(predicted >= cutpoint_nne, 1, 0),
           class_youden = ifelse(predicted >= cutpoint_youden, 1, 0)
           ) |>
    group_by(auc, p0) |>
    slice_sample(n = samples) |>
    arrange(desc(predicted)) |>
    mutate(rank = row_number(),
           classrank5 = ifelse(rank <= 5, 1, 0),
           classrank10 = ifelse(rank <= 10, 1, 0),
           classrank15 = ifelse(rank <= 15, 1, 0),
           classrank20 = ifelse(rank <= 20, 1, 0),
           classrank25 = ifelse(rank <= 25, 1, 0)
           ) |>
    ungroup() |>
    arrange(auc, p0)
}

classifier <- function(classifier, actual, payoffs) {
  case_when(
    classifier == 1 & actual == 1 ~ payoffs[["TP"]],
    classifier == 1 & actual == 0 ~ payoffs[["FP"]],
    classifier == 0 & actual == 1 ~ payoffs[["FN"]],
    classifier == 0 & actual == 0 ~ payoffs[["TN"]]
  )
}

obtain_payoffs <- function(df) {
  payoffs <- sampler()

  df_payoff <- df |>
    mutate(nmb_val_opt = classifier(class_val_opt, actual, payoffs),
           nmb_nne = classifier(class_nne, actual, payoffs),
           nmb_youden = classifier(class_youden, actual, payoffs),
           nmb_classrank5 = classifier(classrank5, actual, payoffs),
           nmb_classrank10 = classifier(classrank10, actual, payoffs),
           nmb_classrank15 = classifier(classrank15, actual, payoffs),
           nmb_classrank20 = classifier(classrank20, actual, payoffs),
           nmb_classrank25 = classifier(classrank25, actual, payoffs)
           ) |>
    group_by(auc, p0) |>
    summarise(across(nmb_val_opt:nmb_classrank25, sum))
  
  return(df_payoff)
}
