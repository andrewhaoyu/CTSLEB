#' Helper function for Super_split. Order qual_tune vector
#'
#' @param x mtx correlation matrix
#' @param y qual_tune vector from helper_qual_tune()
#' @return vector

helper_qual_order <- function(x,y) {
  mtx <- x
  qual_tune <- y
  print("Executing helper_qual_order() ... ")
  prs_tune_order <- order(qual_tune, decreasing=TRUE)

  ix_keep <- prs_tune_order[1]
  for (i in 2:length(prs_tune_order)) {
    if (max(abs(mtx[ix_keep, prs_tune_order[i]])) < 0.98)
      ix_keep <- c(ix_keep, prs_tune_order[i])
  }

  print(paste0(length(ix_keep), ' independent PRS'))
  return(ix_keep)

}
