#' Get the SNP set for estimating the covariance matrix of prior distribution in empirical Bayes method
#'
#' @param x description
#' @param n description
#'
#' @return list
#' @export

PRSTrainValSplit <- function(x,
                             n = 0.50) {
  mat_eb <- x
  n.test <- dim(mat_eb)[1]*n

  super_tune <- as.data.frame(mat_eb[1:n.test,-c(1:2),drop=F])
  super_validate <- as.data.frame(mat_eb[(n.test+1):nrow(mat_eb),-c(1:2),drop=F])

  # drop all the prs columns with pairwise correlation more than 0.98
  # TODO: WE NEED TO MAKE CARET A DEPENDENCY ####
  mtx <- cor(super_tune)

  drop <- findCorrelation(mtx, cutoff=0.98)
  drop <- names(super_tune)[drop]
  super_tune_clean <- super_tune %>%
    select(-all_of(drop))
  super_validate_clean <- super_validate %>%
    select(-all_of(drop))

  return_list <- list("tune" = super_tune_clean,
                      "validate" = super_validate_clean)
  return(return_list)
}
