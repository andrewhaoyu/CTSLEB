#' Align GWAS summary stat for the reference population to the target population.
#' @description
#' Align GWAS summary statistics for the reference population to the target population.
#' The summary statistics data must have similar formats and include the
#' following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect
#' allele. BETA is the regression coefficients for linear regression, log-odds
#' ratio for logistic regression. SE is the standard error for BETA.
#' @param sum_target The GWAS summary statistics for the target population.
#' @param sum_ref The GWAS summary statistics for the reference population.
#' The most commonly used reference population is European.
#' @return data.frame object with ref population GWAS summary statistics
#' aligned with target GWAS summary statistics.
#' @importFrom dplyr left_join mutate select
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data <- "data/"
#' sum_EUR <- fread(paste0(data,"EUR_sumdata.txt"),header=T)
#' sum_AFR <- fread(paste0(data,"AFR_sumdata.txt"),header=T)
#' ref_split_file <- "sum_EUR.txt"
#' target_split_file <- "sum_AFR.txt"
#'
#' AlignSum(sum_target = sum_AFR,
#'          sum_ref = sum_EUR)

AlignSum <- function(sum_target,
                     sum_ref)
  {
  print("executing AlignSum()")
  #match alleles
  sum_ref_select <- sum_ref %>%
    mutate(A1_ref = A1,
           BETA_ref = as.numeric(BETA),
           SE_ref = as.numeric(SE),
           P_ref = as.numeric(P),
           SNP = as.character(SNP)) %>% select(SNP,
           A1_ref,
           BETA_ref,
           SE_ref,
           P_ref)
  sum_target <- sum_target %>%
    mutate(BETA = as.numeric(BETA),
           SE = as.numeric(SE),
           P = as.numeric(P),
           A1 = as.character(A1),
           CHR = as.integer(CHR),
           BP = as.integer(BP),
           SNP = as.character(SNP))

  sum_com <- left_join(sum_target,sum_ref_select,by="SNP")
  idx <- which(sum_com$A1!=sum_com$A1_ref)
  sum_com$BETA_ref[idx] <- -sum_com$BETA_ref[idx]
  sum_com <- sum_com %>% select(-A1_ref)
  print("AlignSum() complete ...")
  return(sum_com)

}

