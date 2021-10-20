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
  # genotype matrix
  geno_matrix <-
    genotype %>%
    dplyr::select(phenotype$SampleNo) %>%
    t()
  # validation (sample size)
  cat(paste0("Sample size of phenotype : ", nrow(phenotype), "\n"))
  cat(paste0("Sample size of genotype : ", nrow(geno_matrix), "\n"))
  if (nrow(phenotype) != nrow(geno_matrix)) {
    stop("Sample size of phenotype is not equal to sample size of genotype")
  }
  # result (algorithm scale)
  result_algo <-
    phenotype %>%
    tidyr::pivot_longer(-SampleNo) %>%
    dplyr::group_by(name) %>%
    tidyr::nest() %>%
    apply(1, as.list) %>%
    purrr::map(GP_cv_repeat,
               genotype_matrix = geno_matrix,
               nfold = nfold,
               nrepeat = nrepeat,
               algorithm = algorithm,
               outputdir = outputdir) %>%
    purrr::flatten()
  # set output list (1st: algorithm name, 2nd: result (algorithm scale))
  output <-
    list(algorithm = algorithm,
         prediction_result = result_algo)
  return(output)
}
