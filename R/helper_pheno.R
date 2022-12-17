#' Helper function for PRS_tune()
#' @description
#' Select file or dataframe source for phenotype y_tun dataframe.
#' Not exported
#' @param pheno_file
#' @param pheno_df
#' @return 'prs_mat' global variable. Matrix of combined PRSs with clumping and
#' pthres as column name
#'
helper_pheno <- function(pheno_file,
                         pheno_df)
  {
  #load the phenotype data for the tuning set
  if (!is.null(pheno_file)) {
    print("phenotype information in ", pheno_file)
    y_tun <- fread(pheno_file)
  } else if (!is.null(pheno_df)) {
    print("phenotype information in a dataframe")
    y_tun <- fread(pheno_df)
  }
  return(y_tun)
}
