#' Generate the combinations of all the tumor characteristics.
#'
#' @param sum_tar The GWAS summary statistics for the target population. The data needs to have following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect allele. BETA is the regression coefficients for linear regression, log-odds ratio for logistic regression. SE is the standard error for BETA.
#' @param sum_other The GWAS summary statistics for the other population. The most commonly used other population is European. The data have similar format as sum_tar.
#'
#' @return GWAS summary statistics for the target population with aligned effects, standard error and p-value for the other population
#' @export
#'
#'
#'

AlignSum = function(sum_tar,sum_other){
  #match alleles
  sum_other_select = sum_other %>%
    mutate(A1_other = A1,
           BETA_other = as.numeric(BETA),
           SE_other = as.numeric(SE),
           P_other = as.numeric(P)，
           SNP = as.numeric(SNP)) %>%
    select(SNP,
           A1_other,
           BETA_other,
           SE_other,
           P_other)
  sum_tar = sum_tar %>%
    mutate(BETA = as.numeric(BETA),
           SE = as.numeric(SE),
           P = as.numeric(P),
           A1 = as.character(A1),
           CHR = as.integer(CHR),
           BP = as.integer(BP),
           SNP = as.character(SNP))

  sum_com <- left_join(sum_tar,sum_other_select,by="SNP")
  idx <- which(sum_com$A1!=sum_com$A1_other)
  sum_com$BETA_other[idx] <- -sum_com$BETA_other[idx]
  sum_com = sum_com %>%
    select(-A1_other)
  return(sum_com)
}
