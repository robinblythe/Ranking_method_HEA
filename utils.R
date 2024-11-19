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