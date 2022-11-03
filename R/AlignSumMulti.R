#' Generate the combinations of all the tumor characteristics.
#'
#' @param sum_tar The GWAS summary statistics for the target population. The data needs to have following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect allele. BETA is the regression coefficients for linear regression, log-odds ratio for logistic regression. SE is the standard error for BETA.
#' @param sum_other_list A list that contains the GWAS summary statistics for multiple populations.
#' @param other_ans_names The names of the other anescestries for GWAS summary statistics
#' @return GWAS summary statistics for the target population with aligned effects, standard error and p-value for the other populations
#' @export
#'
#'
#'

AlignSumMulti <- function(sum_target,sum_other_list,
                         other_ans_names){
  coeff_other_list = list()
  for(i in 1:length(other_ans_names)){
    sum_other_temp  = sum_other_list[[i]]
    sum_com_temp <- AlignSum(sum_target = sum_target,
                             sum_ref = sum_other_temp)
    #a temporary matrix to save the aligned cofficients for the target population
    coeff_other <- sum_com_temp[,c("BETA_ref","SE_ref","P_ref")]
    colnames(coeff_other) <- paste0(c("BETA_","SE_","P_"),other_ans_names[i])
    coeff_other_list[[i]] <- coeff_other
  }
  coeff_mat <- bind_cols(coeff_other_list)

  sum_com <- cbind(sum_target,coeff_mat) %>%
    mutate(P = as.numeric(P),
           BETA = as.numeric(BETA),
           SE = as.numeric(SE))
  return(sum_com)
}
