#' Align GWAS summary stat for the reference population to the target population.
#' @description
#' Align GWAS summary statistics for the reference population to the target population.
#' The summary statistics data must have similar formats and include the
#' following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect
#' allele. BETA is the regression coefficients for linear regression, log-odds
#' ratio for logistic regression. SE is the standard error for BETA.
#' @param sum_tar The GWAS summary statistics for the target population.
#' @param sum_ref The GWAS summary statistics for the reference population.
#' The most commonly used reference population is European.
#' @return data.frame object with ref population GWAS summary statistics
#' aligned with target GWAS summary statistics. The resulting global variable is
#' named "sum_com".
#' @export
#' @examples
#' data.dir <- "data/"
#' temp.dir = "test/temp.dir/"
#' sum_EUR <- fread(paste0(data.dir,"EUR_sumdata.txt"),header=T)
#' sum_AFR <- fread(paste0(data.dir,"AFR_sumdata.txt"),header=T)
#' ref_split_file <- "sum_EUR.txt"
#' target_split_file <- "sum_AFR.txt"
#'
#' AlignSum(sum_target = sum_AFR,
#'          sum_ref = sum_EUR,
#'          results_dir = results.dir,
#'          ref_split_file = ref_split_file,
#'          target_split_file = target_split_file,
#'          SplitSum = TRUE)

AlignSum <- function(sum_target,
                     sum_ref)
  {

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

  sum.com <- left_join(sum_target,sum_ref_select,by="SNP")
  idx <- which(sum.com$A1!=sum.com$A1_ref)
  sum.com$BETA_ref[idx] <- -sum.com$BETA_ref[idx]
  sum.com <- sum.com %>% select(-A1_ref)

  # if (SplitSum) {
  #   #assign("sum_com", sum.com, envir = .GlobalEnv)
  #
  #   list <- helper_SplitSum(sum_com = sum.com,
  #                           results_dir = results_dir)
  #   list <- c(list, sum_com=sum.com)
  #   return(list)
  # } else {
  #   print(paste0("SplitSum() was not performed"))
  #   return(sum.com)
  # }
}

