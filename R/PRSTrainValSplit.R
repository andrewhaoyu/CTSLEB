#' Get the SNP set for estimating the covariance matrix of prior distribution in empirical Bayes method
#'
#' @param x matrix of PRSs based on EB coefficients produced by PRSscoreEBayes()
#' or CalculateEBEffectSize()
#' @param n fraction to use for validation
#'
#' @return list
#' @export

PRSTrainValSplit <- function(x,
                             n = 0.5) {
  print("Executing PRSTrainValSplit() ... ")
  mat_eb <- x
  n.test <- dim(mat_eb)[1]*n

  super_tune <- as.data.frame(mat_eb[1:n.test,-c(1:2),drop=F])
  super_validate <- as.data.frame(mat_eb[(n.test+1):nrow(mat_eb),-c(1:2),drop=F])

  # drop all the prs columns with pairwise correlation more than 0.98
  # TODO: WE NEED TO MAKE CARET A DEPENDENCY ####
  print("Executing correlation ... ")
  mtx <- cor(super_tune)
  assign("mtx", mtx, envir = .GlobalEnv)
  drop <- findCorrelation(mtx, cutoff=0.98)
  drop <- names(super_tune)[drop]
  super_tune_clean <- super_tune %>%
    select(-all_of(drop))
  super_validate_clean <- super_validate %>%
    select(-all_of(drop))

  return_list <- list("tune" = super_tune_clean,
                      "validate" = super_validate_clean)
  print("tune object created ... ")
  print("validate object created ... ")
  return(return_list)
}

