GPwithGWAS_core <- function(list,
                            nSNPs,
                            nfold,
                            nrepeat,
                            algorithm,
                            method,
                            outputdir) {

  phenotype_name <- list$name
  phenotype_dir <- paste0(outputdir, "/", algorithm, "/", method, "/", phenotype_name)
  fs::dir_create(phenotype_dir)

  y <- list$y

  output <- list()

  for (i in seq_len(length(nSNPs))) {
    n_snp <- nSNPs[i]
    cat(insight::print_color(paste0("nSNPs      : ", n_snp, "\n"), "green"))

    x <-
      list$geno %>%
      dplyr::slice(1:n_snp) %>%
      dplyr::select(y$key) %>%
      t()

    result_snp <-
      GP_cv_repeat(
        phenotype_name = phenotype_name,
        phenotype_data = y,
        genotype_matrix = x,
        nfold = nfold,
        nrepeat = nrepeat,
        algorithm = algorithm,
        outputdir = outputdir,
        output = FALSE) %>%
      `names<-` (NULL) %>%
      as.data.frame()

    utils::write.csv(result_snp,
                     file = paste0(phenotype_dir, "/GP_", n_snp, "SNPs_", nfold, "-fold_", nrepeat, "-repeat_Prediction.csv"))

    output <- append(output, list(result_snp))
    names(output)[i] <- paste0(n_snp, "SNPs")
  }
  output <- list(output)
  names(output) <- phenotype_name

  return(output)
}
