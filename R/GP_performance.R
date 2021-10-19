#' Calculation of R and RMSE
#'
#' @importFrom dplyr %>%
#'
#' @param prediction_df dataframe of prediction result
#'
GP_evaluation <- function(prediction_df) {

  predicted_y <-
    prediction_df %>%
    dplyr::select(starts_with("Rep"))
  observed_y <-
    prediction_df %>%
    dplyr::select(-starts_with("Rep")) %>%
    as.matrix()
  # get phenotype name
  this_pheno <- colnames(observed_y)
  # R
  r <- apply(predicted_y, 2, function(x) {stats::cor(x, observed_y)})
  # R2
  r2 <- apply(predicted_y, 2, function(x) {MLmetrics::R2_Score(x, observed_y)})
  # RMSE
  rmse <- apply(predicted_y, 2, function(x) {MLmetrics::RMSE(x, observed_y)})
  # set output list (R, R2 & RMSE)
  output <- list(
    R = data.frame(r_df = r) %>% dplyr::rename(!!this_pheno := r_df),
    R2 = data.frame(r2_df = r2) %>% dplyr::rename(!!this_pheno := r2_df),
    RMSE = data.frame(rmse_df = rmse) %>% dplyr::rename(!!this_pheno := rmse_df)
  )
  return(output)
}

#' Evaluating the result of TeaWAS::GP_cv
#'
#' @importFrom dplyr %>%
#'
#' @param prediction_list List of return values for TeaWAS::GP_cv
#'
#' @export
#'
GP_performance <- function(prediction_list) {
  # get algorithm name
  this_algorithm <- prediction_list$algorithm
  # evaluation
  performance_list <-
    prediction_list$prediction_result %>%
    purrr::map(GP_evaluation)
  # set output list (R, R2 & RMSE)
  output <-
    list(
      R =
        performance_list %>%
        purrr::map(~ {use_series(., R)}) %>%
        purrr::reduce(dplyr::bind_cols),
      R2 =
        performance_list %>%
        purrr::map(~ {use_series(., R2)}) %>%
        purrr::reduce(dplyr::bind_cols),
      RMSE =
        performance_list %>%
        purrr::map(~ {use_series(., RMSE)}) %>%
        purrr::reduce(dplyr::bind_cols)
    ) %>%
    purrr::map(~ {
        .x %>%
        dplyr::mutate(
          n = as.numeric(stringr::str_extract_all(rownames(.), "[0-9]")),
          algorithm = this_algorithm)
    })
  return(output)
}
