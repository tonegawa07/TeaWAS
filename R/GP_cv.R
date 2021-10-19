#' Genomic Prediction by k-fold cross validation
#'
#' @importFrom dplyr %>%
#'
#' @param phenotype dataframe of phenotype
#' @param genotype matrix of SNPs data
#' @param nfold Number of divisions for cross-validation
#' @param nrepeat Number of prediction iterations
#' @param algorithm Algorithm to be used
#' @param outputdir Directory to output results
#'
#' @export
#'
GP_cv <- function(phenotype, genotype, nfold, nrepeat, algorithm, outputdir) {
  # get algorithn name
  this_algorithm <- algorithm
  # make algorithm name dir
  dir.create(paste0(outputdir, "/", this_algorithm))
  # set output list (1st: algorithm name)
  output <- list(algorithm = this_algorithm)
  # result (algorithm scale)
  result_algo <-  list()

  # for loop (phenotype)
  for (n_pheno in seq_len(ncol(phenotype))) {
    # get phenotype name
    this_pheno <- colnames(phenotype)[n_pheno]
    # print information
    print("", quote = F)
    print(paste0("algorithm : ", this_algorithm), quote = F)
    print(paste0("phenotype : No.", n_pheno, " ", this_pheno), quote = F)
    # make phenotype name dir
    dir.create(paste0(outputdir, "/", this_algorithm, "/", this_pheno))
    # get objective variable
    y <-
      phenotype %>%
      dplyr::select(this_pheno)
    # make dataset (bind y & x)
    dataset <- cbind(y, genotype)
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
        purrr::map(GP_cv_core, phenotype_name = this_pheno, algorithm = this_algorithm) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::arrange(rownames(.)) %>%
        dplyr::rename(!!paste0("Rep_", n_rep) := pred)
      # result (phenotype scale)
      result_pheno <-
        result_pheno %>%
        dplyr::bind_cols(result_rep)
    } # for loop (repeat)
    # output result (phenotype scale)
    utils::write.csv(result_pheno,
              file = paste0(outputdir, "/", this_algorithm, "/", this_pheno, "/GP_", nfold, "-fold_", nrepeat, "-repeat_Prediction.csv")
              )
    # result (algorithm scale)
    result_algo <- c(result_algo, pheno_name = list(result_pheno))
    names(result_algo)[n_pheno] <- this_pheno
  } # for loop (phenotype)
  # set output list (2nd: result (algorithm scale))
  output <- c(output, prediction_result = list(result_algo))
  return(output)
}
