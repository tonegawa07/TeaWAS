GWAS_multiple <- function(phenotype,
                          genotype,
                          MAF = 0.05,
                          K = FALSE,
                          Q_npc = 0,
                          outputdir) {

  # phenotype dataframe
  phenotype <-
    phenotype %>%
    stats::na.omit()

  genotype <-
    genotype %>%
    dplyr::select(id = id, chr = chr, pos = pos, rownames(phenotype))
  # samplesize of geno
  ngeno <- genotype %>% dplyr::select(rownames(phenotype)) %>% ncol()

  # validation (sample size)
  cat(paste0("Sample size of phenotype : ", nrow(phenotype), "\n"))
  cat(paste0("Sample size of genotype : ", ngeno, "\n"))
  if (nrow(phenotype) != ngeno) {
    stop("Sample size of phenotype is not equal to sample size of genotype")
  }

  if (K) {
    this_K <-
      genotype %>%
      dplyr::select(rownames(phenotype)) %>%
      t() %>%
      rrBLUP::A.mat(shrink = TRUE)
  } else {
    this_K <- NULL
  }

  params <-
    list(
      MAF.  = MAF,
      K.    = this_K,
      n.PC. = Q_npc
    )

  if (is.null(params$K.) && params$n.PC.==0) {
    model <- "GLM"
  } else if (is.null(params$K.)) {
    model <- paste0("MLM_Q_pc", params$n.PC.)
  } else if (params$n.PC.==0) {
    model <- "MLM_K"
  } else {
    model <- paste0("MLM_QK_pc", params$n.PC.)
  }
  cat(insight::print_color(paste0("Model : ", model, "\n"), "green"))

  modeldir <- paste0(outputdir, "/", model)

  output <-
    phenotype %>%
    dplyr::mutate(key = rownames(.)) %>%
    tidyr::pivot_longer(-key) %>%
    dplyr::group_by(name) %>%
    tidyr::nest() %>%
    apply(1, as.list) %>%
    purrr::map(GWAS_core,
               genotype = genotype,
               params = params,
               outputdir = modeldir) %>%
    purrr::flatten()
}
