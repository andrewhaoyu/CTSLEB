#' Helper function for PreparePlinkFile()
#' @description
#' Calls CreateQRange() using either the default pthres or user input values
#' Not exported
#' @param x pthres
#' @return CreateQRange() object
#'
helper_return_list <- function(x,
                               score_file,
                               p_value_file,
                               unique_infor,
                               q_range)
  {
  if (x) {
    result <- list(score_file,
                   p_value_file,
                   unique_infor,
                   q_range)
    return(result)
  } else {
    assign("score_file", score_file, envir = .GlobalEnv)
    assign("p_value_file", p_value_file, envir = .GlobalEnv)
    assign("unique_infor", unique_infor, envir = .GlobalEnv)
    assign("q_range", q_range, envir = .GlobalEnv)
  }
}
