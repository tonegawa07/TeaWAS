GPwithGWAS_evaluation <- function(list) {
  output <-
    list %>%
    purrr::map(GP_evaluation) %>%
    purrr::map(~ {magrittr::use_series(., R)}) %>%
    purrr::transpose() %>%
    `names<-` (NULL) %>%
    as.data.frame() %>%
    dplyr::rename_at(dplyr::vars(dplyr::starts_with("X")),
                     ~ stringr::str_replace(., pattern = "X", replacement = ""))
}

#' @export
GPwithGWAS_performance <- function(list) {

  output <-
    list %>%
    magrittr::use_series(prediction_result) %>%
    purrr::map(GPwithGWAS_evaluation) %>%
    as.matrix() %>%
    as.data.frame() %>%
    dplyr::mutate(name = rownames(.)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(data = purrr::map2(.x = V1, .y = name,
                                     ~dplyr::mutate(.x,
                                                    phenotype = .y,
                                                    n = 1:nrow(.x)))) %>%
    dplyr::select(data) %>%
    as.list() %>%
    magrittr::use_series(data) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    tidyr::pivot_longer(cols = -c(phenotype, n), names_to = "nSNPs") %>%
    tidyr::pivot_wider(names_from = phenotype, values_from = value) %>%
    dplyr::mutate(
      nSNPs = stringr::str_replace_all(nSNPs, c("SNPs" = "")) %>% as.integer(),
      algorithn = list$algorithm
      ) %>%
    dplyr::arrange(nSNPs)

  return(output)

}
