#' Prepare the files for PLINK2 to calculate PRSs
#' @description
#' Calls CreateQRange() using either the default pthres or user input values
#' @param snp_list the snp_list result from the two-dimensional clumping
#' @param sum_com the sum_com result from the AlignSum function.
#' @param pthres vector of p-value thresholds. Default
#' c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#' @param return_list Return as a list instead of creating a global variable
#' @return the score_file, p_value_file and p_other_file for the two-dimensional thresholding
#' @export
#'
#' @examples
PreparePlinkFile <- function(snp_list,
                             sum_com,
                             pthres = as.null(),
                             output,
                             return_list = FALSE)
  {
  #create unique SNP list by combind LD clumping results under different parameters
  unique_id <- unique(rbindlist(snp_list,use.name =FALSE))
  names(unique_id) <- "SNP"
  #align the regression coefficients for these SNPs from the sum stat
  unique_infor <- left_join(unique_id,sum_com,by="SNP")

  #create a coefficient matrix
  #the first column contains the unique SNPs after clumping results under all combinations of r2-cutoff and window_size
  #the second column is the effect allele
  #the third to the last columns contains the regression coefficients of the target population for SNPs after LD-clumping under a specific combination of r2-cutoff and base_window_size
  #the coefficients is put as 0 if a SNP doesn't exist in the clumping results under a specific combination of r2-cutoff and base_window_size
  n_col <- length(snp_list)
  n_row <- nrow(unique_infor)
  beta_mat <- matrix(unique_infor$BETA,nrow =n_row,ncol =n_col)
  names <- rep("c",n_col)
  temp <- 1
  for(ldx in seq(n_col)){
    LD <- snp_list[[ldx]]
    names(LD) <- "SNP"
    idx <- which(unique_infor$SNP%in%LD$SNP==F)
    beta_mat[idx,ldx] <- 0
    names[ldx] <- names(snp_list[[ldx]])
  }
  colnames(beta_mat) <- names
  score_file <- data.frame(SNP = unique_id,A1 = unique_infor$A1,beta_mat)
  print("score_file complete")

  p_value_file <- data.frame(SNP = unique_id,P = unique_infor$P)
  helper_p_value_file(x = p_value_file,
                      pthres = pthres,
                      output = output)
  print("p_value_file complete")

  q_range <- helper_CreateQRange(pthres)
  print("q_range complete")

  helper_return_list(x = return_list,
                     output = output,
                     score_file,
                     p_value_file,
                     unique_infor,
                     q_range)
#  if (is.null(pthres)) {
#    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,
#                5E-02,5E-01,1.0)
#    print(paste0("pthres is NULL...using default values "))
#    q_range <- CreateQRange(pthres)
#  }

  # if (return_list) {
  #   result <- list(score_file,
  #                  p_value_file,
  #                  unique_infor)
  # } else {
  #   assign("score_file", score_file, envir = .GlobalEnv)
  #   assign("p_value_file", p_value_file, envir = .GlobalEnv)
  #   assign("unique_infor", unique_infor, envir = .GlobalEnv)
  # }
  #
  # assign("q_range", q_range, envir = .GlobalEnv)
  # return(result)

}
