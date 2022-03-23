#' Generate the combinations of all the tumor characteristics.
#'
#' @param sum_tar The GWAS summary statistics for the target population. The data needs to have following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect allele. BETA is the regression coefficients for linear regression, log-odds ratio for logistic regression. SE is the standard error for BETA.
#' @param sum_other The GWAS summary statistics for the other population. The most commonly used other population is European. The data have similar format as sum_tar.
#'
#' @return GWAS summary statistics for the target population with aligned effects, standard error and p-value for the other population
#' @export
#'
#' @examples
#'
#'#this is a simulated breast cancer example
#'#there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
#' data[1:5,]

#'#four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
#'#the phenotype file
#'y <- data[,1:5]
#'#generate the combinations of all the subtypes
#'#by default, we remove all the subtypes with less than 10 cases
#'z.standard <- GenerateZstandard(y)

AlignSum <- function(sum_tar,sum_other){
  #match alleles
  sum_other_select = sum_other %>%
    rename(A1_other = A1,
           BETA_other = BETA,
           SE_other = SE,
           P_other = P) %>%
    select(SNP,
           A1_other,
           BETA_other,
           SE_other,
           P_other)
  sum_com <- left_join(sum_tar,sum_other_select,by="SNP")
  idx <- which(sum_com$A1!=sum_com$A1_other)
  sum_com$BETA_other[idx] <- -sum_com$BETA_other[idx]
  sum_com = sum_com %>%
    select(-A1_other)
  return(sum_com)
}
