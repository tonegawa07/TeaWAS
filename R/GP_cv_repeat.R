#' Genomic Prediction per phenotype
#'
#' @importFrom dplyr %>%
#'
#' @param phenotype_name name of each phenotype
#' @param phenotype_data dataframe of each phenotype
#' @param genotype_matrix Transposed genotype matrix with equal number of rows
#' @param nfold Number of divisions for cross-validation
#' @param nrepeat Number of prediction iterations
#' @param algorithm Algorithm to be used
#' @param outputdir Directory to output results
#'
GP_cv_repeat <- function(phenotype_name,
                         phenotype_data,
                         genotype_matrix,
                         nfold,
                         nrepeat,
                         algorithm,
                         outputdir) {
  # print information
  cat(insight::print_color(paste0("Algorithm : ", algorithm, "\n"), "green"))
  cat(insight::print_color(paste0("Phenotype : ", phenotype_name, "\n"), "green"))
  # make phenotype name dir
  dir.create(paste0(outputdir, "/", algorithm, "/", phenotype_name))
  # get objective variable
  y <-
    phenotype_data %>%
    as.data.frame() %>%
    `rownames<-`(.$key) %>%
    dplyr::select(!!phenotype_name := value)
  # make dataset (bind y & x)
  dataset <- cbind(y, genotype_matrix)
  # result (phenotype scale)
  result_pheno <- y
  # set progress bar params
  pb <- utils::txtProgressBar(min = 0, max = nrepeat, style = 3)

  # for loop (repeat)
  for (n_rep in 1:nrepeat) {
    # progress bar
    utils::setTxtProgressBar(pb, n_rep)
    # result (repeat scale)
    result_rep <-
      dataset %>%
      rsample::vfold_cv(v = nfold) %>%
      magrittr::use_series(splits) %>%
      purrr::map(GP_cv_core,
                 phenotype_name = phenotype_name,
                 algorithm = algorithm) %>%
      purrr::reduce(dplyr::bind_rows) %>%
      dplyr::arrange(rownames(.)) %>%
      dplyr::rename(!!paste0("Rep_", n_rep) := pred)
    # result (phenotype scale)
    result_pheno <-
      result_pheno %>%
      dplyr::bind_cols(result_rep)
  } # for loop (repeat)
  # New line in console
  cat("\n")
  # output result (phenotype scale)
  utils::write.csv(result_pheno,
                   file = paste0(outputdir, "/", algorithm, "/", phenotype_name, "/GP_", nfold, "-fold_", nrepeat, "-repeat_Prediction.csv"))
  # set output list (result (phenotype scale))
  output <- list(result_pheno)
  names(output) <- phenotype_name

  return(output)
}
