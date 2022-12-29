#' Prepare the files for PLINK2 to calculate PRSs
#' @description
#' Calls helper_return_list() and helper_CreateQRange()
#' @param snp_list the snp_list result from the two-dimensional clumping
#' @param sum_com the sum_com result from the AlignSum function.
#' @param pthres vector of p-value thresholds. Default
#' c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#' @param results_dir
#' @return Creates either the global variables 'scores', 'p_values',
#' 'q_range' and 'unique_infor' or the variable 'plink_list' which contains the
#' four dataframes in a list
#' @export
#' @examples
PreparePlinkFile <- function(snp_list = snp_list,
                             sum_com,
                             pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                             results_dir,
                             params_farm=as.null())
  {
  print("executing PreparePlinkFile()... ")
  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }

  #create unique SNP list by combind LD clumping results under different parameters
  unique_id <- unique(rbindlist(snp_list,use.name =FALSE))
  names(unique_id) <- "SNP"

  #align the regression coefficients for these SNPs from the sum stat

  unique_infor <- left_join(unique_id,sum_com,by="SNP")
  print("unique_infor dataframe complete")

  #create a coefficient matrix

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

  scores <- data.frame(SNP = unique_id,A1 = unique_infor$A1,beta_mat)
  print("scores dataframe complete")
  p_values <- data.frame(SNP = unique_id,P = unique_infor$P)
  print("p_values dataframe complete")
  q_range <- CreateQRange(pthres)

  names <- c("scores",
             "p_values",
             "unique_infor",
             "q_range")
  values <- list(scores,
                 p_values,
                 unique_infor,
                 q_range)
  list <- setNames(values, names)
  print("PreparePlinkFile() complete ...")

  return(list)

}
