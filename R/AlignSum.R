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
#' @param SplitSum Execute the SplitSum() and WriteSplitTables() function.
#' Default = TRUE.
#' @param ref_split_file If SplitSum = TRUE (Default) provide a file path and prefix name
#' for the reference snps produced from the SplitSum() function.  The default
#' name "sum" is passed to WriteSplitTables() which creates the file
#' "sum_ref.txt" in the current working directory. This file is the input for the
#' plink --clump flag
#' @param target_split_file If SplitSum = TRUE (Default) provide a file path and prefix name
#' for the target snps produced from the SplitSum() function. The default
#' name "sum" is passed to WriteSplitTables() which creates the file
#' "sum_target.txt" in the current working directory. This file is the input for
#' the plink --clump flag
#' @param results.dr
#' @return data.frame object with ref population GWAS summary statistics
#' aligned with target GWAS summary statistics. The resulting global variable is
#' named "sum_com". If SplitSum = TRUE, the SplitSum() and WriteSplitTables()
#' functions will be executed and the following global variables will be created;
#' "sum_ref" for the reference snps and "sum_target" for the target snps created
#' by SplitSum(). The data.frames are also written to files named "sum_ref.txt"
#' and "sum_target.txt" respectively
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
                     sum_ref,
                     results_dir = results.dir,
                     ref_split_file="sum_ref.txt",
                     target_split_file="sum_target.txt",
                     SplitSum=TRUE
                     )
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

  sum_com <- left_join(sum_target,sum_ref_select,by="SNP")
  idx <- which(sum_com$A1!=sum_com$A1_ref)
  sum_com$BETA_ref[idx] <- -sum_com$BETA_ref[idx]
  sum_com <- sum_com %>% select(-A1_ref)

  assign("sum_com", sum_com, envir = .GlobalEnv)

  if (SplitSum) {
    helper_SplitSum(sum_com = sum_com,
                    ref_split_file = ref_split_file,
                    target_split_file = target_split_file)
  } else {
    print(paste0("SplitSum() was not performed"))
  }
}

