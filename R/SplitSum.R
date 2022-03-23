#' Generate the combinations of all the tumor characteristics.
#'
#' @param sum_com The GWAS summary statistics for the target population. The data needs to have following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect allele. BETA is the regression coefficients for linear regression, log-odds ratio for logistic regression. SE is the standard error for BETA.
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

SplitSum <- function(sum_com){
  sum_com_select = sum_com %>%
    mutate(split_ind =
             ifelse(
               (P<P_other)|is.na(P_other),T,F)
           )%>%
    select(SNP,P,P_other,split_ind)

  sum_com_select_other_ref = sum_com_select %>%
    filter(split_ind==F) %>%
    select(SNP,P_other) %>%
    rename(P = P_other)

  sum_com_select_tar_ref = sum_com_select %>%
    filter(split_ind==T) %>%
    select(SNP,P)

  sum_list = list(sum_com_select_other_ref,
                  sum_com_select_tar_ref)
  return(sum_list)

}
