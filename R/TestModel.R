#' Get the SNP set for estimating the covariance matrix of prior distribution in empirical Bayes method
#'
#' @param x description
#' @param y description
#' @param y validate
#' @param ml_library description
#'
#' @return lm model
#' @export

TestModel <- function(x,y,
                      validate,
                      ml_library = c(
                        "SL.glmnet",
                        "SL.ridge",
                        "SL.nnet" )) {
  this_y <- y
  this_x <- x

  sl <- SuperLearner(Y = this_y, X = this_x, family = gaussian(),
                     SL.library = ml_library)

  #predict the outcome using the independent validation dataset

  y_pred <- predict(sl, validate, onlySL = TRUE)

  #evaluate the CT-SLEB prs performance on the validation

  mod <- lm(validate~y_pred[[1]])
  return(mod)
}
