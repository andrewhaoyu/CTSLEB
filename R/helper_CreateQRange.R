#' Helper function for PreparePlinkFile()
#' @description
#' Calls CreateQRange() with either the default pthres or user input values
#' Not exported
#' @param x pthres
#' @return CreateQRange() object
#'
helper_CreateQRange <- function(x)
  {
  if (is.null(x)) {
    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,
                5E-02,5E-01,1.0)
    print(paste0("pthres is NULL...using default values "))
  } else {
    pthres <- x
  }
  print(x)
  print(pthres)
  q_range <- CreateQRange(pthres)
  return (q_range)
}
