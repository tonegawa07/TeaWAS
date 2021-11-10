GWAS_core <- function(phenotype_list, genotype, params, outputdir) {
  # get phenotype name
  this_pheno <- phenotype_list$name

  pheno. <-
    phenotype_list$data %>%
    as.data.frame() %>%
    `rownames<-`(.$key) %>%
    dplyr::select(gid = key, !!this_pheno := value)

  geno. <- genotype

  result_pheno <-
    rrBLUP::GWAS(pheno = pheno.,
                 geno = geno.,
                 K = params$K.,
                 n.PC = params$n.PC.,
                 min.MAF = params$MAF.,
                 plot = FALSE) %>%
    dplyr::rename(y := !!this_pheno) %>%
    dplyr::filter(y != 0) %>%
    dplyr::mutate(P = 10^(-y)) %>%
    dplyr::select(SNP = id, CHR = chr, BP = pos, P)

  fs::dir_create(c(paste0(outputdir, "/table"), paste0(outputdir, "/fig")))
  utils::write.csv(result_pheno,
                   file = paste0(outputdir, "/table/", this_pheno,  "_gwas_result.csv"),
                   row.names = FALSE)

  #draw manhattan plot
  png(paste0(outputdir, "/fig/", this_pheno, "_manhattan.png"), width = 1200, height = 732)
  qqman::manhattan(result_pheno,
                   cex.lab = 2,
                   cex.axis = 2,
                   cex.main = 2,
                   cex = 2,
                   suggestiveline = FALSE,
                   genomewideline = FALSE,
                   main = this_pheno,
                   col = c("dodgerblue4", "gray66"))
  dev.off()

  ##draw Q-Q plot
  png(paste0(outputdir, "/fig/", this_pheno, "_QQplot.png"), width = 300, height = 300)
  qqman::qq(result_pheno$P,
            main = this_pheno,
            col = "black",
            cex.lab = 1,
            cex.axis = 1,
            cex.main = 1,
            cex = 1)
  dev.off()

  # set output list (result (phenotype scale))
  output <- list(result_pheno)
  names(output) <- this_pheno

  return(output)
}
