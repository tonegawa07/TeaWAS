GPwithGWAS <- function(phenotype, genotype, GWAS_result, nSNPs, nfold, nrepeat, algorithm, method, outputdir) {

  phenotype. <-
    phenotype %>%
    dplyr::mutate(key = rownames(.)) %>%
    tidyr::pivot_longer(-key) %>%
    dplyr::group_by(name) %>%
    tidyr::nest() %>%
    dplyr::rename(y = data)

  genotype. <-
    GWAS_result %>%
    tidyr::pivot_longer(cols = -c(1:3), values_to = "P") %>%
    dplyr::group_by(name) %>%
    tidyr::nest() %>%
    dplyr::mutate(selected_SNPs = purrr::map(data,
                                             ~dplyr::slice(
                                               dplyr::arrange(.x, P),
                                               1:max(nSNPs)))
                  ) %>%
    dplyr::select(-data) %>%
    dplyr::mutate(geno = purrr::map(selected_SNPs,
                                    ~dplyr::inner_join(.x, genotype)))

  dataset <-
    genotype. %>%
    dplyr::inner_join(phenotype.) %>%
    dplyr::select(name, geno, y) %>%
    apply(1, as.list)

  result_algo <-
    dataset %>%
    purrr::map(GPwithGWAS_core,
               nSNPs = nSNPs,
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
