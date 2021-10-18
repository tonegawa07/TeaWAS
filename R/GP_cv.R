#' Genomic Prediction by k-fold cross validation
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
GP_cv <- function(phenotype, genotype, nfold, nrepeat, algorithm, outputdir)
{
  # get algorithn name
  this_algorithm <- algorithm
  # make algorithm name dir
  dir.create(paste0(outputdir, "/", this_algorithm))
  # set output list (1st: algorithm name)
  output <- list(algorithm = this_algorithm)
  # result (algorithm scale)
  result_algo <-  list()

  # for loop (phenotype)
  for (n_pheno in 1:ncol(phenotype)) {

    this_pheno <- colnames(pheno)[n_pheno]
    print("", quote = F)
    print(paste0("algorithm : ", this_algorithm), quote = F)
    print(paste0("phenotype : No.", n_pheno, " ", this_pheno), quote = F)

    dir.create(paste0(outputdir, "/", this_algorithm, "/", this_pheno))

    y <-
      phenotype %>%
      dplyr::select(this_pheno)

    dataset <- cbind(y, genotype)

    result_pheno <- y

    pb <- txtProgressBar(min = 0, max = nrepeat, style = 3)

    # for loop (repeat)
    for (n_rep in 1:nrepeat) {

      setTxtProgressBar(pb, n_rep)

      result_rep <-
        dataset %>%
        rsample::vfold_cv(v = nfold) %>%
        magrittr::use_series(splits) %>%
        purrr::map(GP, phenotype_name = this_pheno, algorithm = this_algorithm) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::arrange(rownames(.)) %>%
        dplyr::rename(!!paste0("Rep_", n_rep) := pred)

      result_pheno <-
        result_pheno %>%
        dplyr::bind_cols(result_rep)
    } # for loop (repeat)

    write.csv(result_pheno,
              file = paste0(outputdir, "/", this_algorithm, "/", this_pheno, "/GP_", nfold, "-fold_", nrepeat, "-repeat_Prediction.csv")
              )

    result_algo <- c(result_algo, pheno_name = list(result_pheno))
    names(result_algo)[n_pheno] <- this_pheno
  } # for loop (phenotype)

  output <- c(output, prediction_result = list(result_algo))
  return(output)
}
