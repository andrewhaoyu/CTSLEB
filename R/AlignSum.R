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
#' @param SplitSum Execute the SplitSum() function. Default = TRUE.
#' @return data.frame object with ref population GWAS summary statistics
#' aligned with target GWAS summary statistics. The resulting global variable is
#' named "sum_com". If SplitSum = TRUE, the SplitSum() and WriteSplitTables()
#' functions will be executed and the following global variables will be created;
#' "sum_ref" for the reference snps and "sum_target" for the target snps created
#' by SplitSum().
#' @export
#' @examples
#' data.dir <- "data/"
#' sum_EUR <- fread(paste0(data.dir,"EUR_sumdata.txt"),header=T)
#' sum_AFR <- fread(paste0(data.dir,"AFR_sumdata.txt"),header=T)
#' AlignSum(sum_tar = sum_AFR, sum_ref = sum_EUR, SplitSum=TRUE)

AlignSum <- function(sum_target,sum_ref,SplitSum=TRUE){
  #match alleles
  sum_ref_select <- sum_ref %>%
    mutate(A1_ref = A1,
           BETA_ref = as.numeric(BETA),
           SE_ref = as.numeric(SE),
           P_ref = as.numeric(P),
           SNP = as.character(SNP)) %>%
    select(SNP,
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
  sum_com <- sum_com %>%
    select(-A1_ref)

  if (SplitSum) {
    assign("sum_com", sum_com, envir = .GlobalEnv)
    split_list <- SplitSum(sum_com)
    WriteSplitTables(x = split_list)
  } else {
    print(paste0("SplitSum() was not performed"))
    assign("sum_com", sum_com, envir = .GlobalEnv)
  }
}

