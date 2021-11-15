#' Genomic Prediction by k-fold cross validation
#'
#' @importFrom dplyr %>%
#'
#' @param phenotype dataframe of phenotype. rowname is a key number
#' @param genotype matrix of SNPs data. rowname is a key number
#' @param nfold Number of divisions for cross-validation
#' @param nrepeat Number of prediction iterations
#' @param algorithm Algorithm to be used
#' @param outputdir Directory to output results
#'
#' @export
#'
GP_cv <- function(phenotype,
                  genotype,
                  nfold,
                  nrepeat,
                  algorithm,
                  outputdir) {
  # make algorithm name dir
  dir.create(paste0(outputdir, "/", algorithm))
  # phenotype dataframe
  phenotype <- phenotype %>% stats::na.omit()

  # validation (sample size)
  cat(paste0("Sample size of phenotype : ", nrow(phenotype), "\n"))
  cat(paste0("Sample size of genotype : ", nrow(genotype), "\n"))
  if (nrow(phenotype) != nrow(genotype)) {
    stop("Sample size of phenotype is not equal to sample size of genotype")
  }
  # result (algorithm scale)
  result_algo <-
    phenotype %>%
    dplyr::mutate(key = rownames(.)) %>%
    tidyr::pivot_longer(-key) %>%
    dplyr::group_by(name) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      GP = purrr::map2(name, data,
                       ~GP_cv_repeat(
                          phenotype_name = .x,
                          phenotype_data = .y,
                          genotype_matrix = genotype,
                          nfold = nfold,
                          nrepeat = nrepeat,
                          algorithm = algorithm,
                          outputdir = outputdir))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(GP) %>%
    purrr::flatten() %>%
    purrr::flatten()

  # set output list (1st: algorithm name, 2nd: result (algorithm scale))
  output <-
    list(algorithm = algorithm,
         prediction_result = result_algo)
  return(output)
}
