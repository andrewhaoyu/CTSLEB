#' Generate the combinations of all the tumor characteristics.
#'
#' @param sum_com The GWAS summary statistics for the target population. The data needs to have following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect allele. BETA is the regression coefficients for linear regression, log-odds ratio for logistic regression. SE is the standard error for BETA.
#'
#' @return GWAS summary statistics for the target population with aligned effects, standard error and p-value for the other population
#' @export
#'
#' @examples

SplitSum <- function(sum_com){
  sum_com_select <- sum_com %>%
    mutate(split_ind =
             ifelse(
               (P<P_ref)|is.na(P_ref),T,F)
    )%>%
    select(SNP,P,P_ref,split_ind)

  sum_com_select_other_ref <- sum_com_select %>%
    filter(split_ind==F) %>%
    select(SNP,P_ref) %>%
    rename(P = P_ref)

  sum_com_select_tar_ref <- sum_com_select %>%
    filter(split_ind==T) %>%
    select(SNP,P)

  sum_list <- list(sum_com_select_other_ref,
                   sum_com_select_tar_ref)

  return(sum_list)

}
