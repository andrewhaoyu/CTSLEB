#' Helper function for PreparePlinkFile()
#' @description
#' Return either a list containing the score_file, p_values, unique_infor and
#' q_range or the global variables for the same objects and write the files.
#' Not exported
#' @param x return_list
#' @param results_dir
#' @param scores
#' @param p_values
#' @param unique_infor
#' @param q_range
#' @return CreateQRange() object
#'
helper_return_list <- function(x,
                               results_dir,
                               scores,
                               p_values,
                               unique_infor,
                               q_range)
  {

  if (x) {
    print("return plink_file")
    result <- list(scores = scores,
                   p_values = p_values,
                   unique_infor = unique_infor,
                   q_range = q_range)
    assign("plink_files", result, envir = .GlobalEnv)
  } else {
    temp.dir <- paste0(results_dir,"temp/")
    scores_file_out <- paste0(temp.dir, "scores_file")
    print(paste0("printing scores_file: ", scores_file_out))
    write.table(scores,
               file = scores_file_out,
               row.names = F,
               col.names = F,
               quote=F)
    q_range_out <- paste0(temp.dir, "q_range_file")
    print(paste0("printing q_range_file: ", q_range_out))
    write.table(q_range,
               file = q_range_out,
               row.names = F,
               col.names = F,
               quote=F)
    n_col <- ncol(scores)
    assign("scores", scores, envir = .GlobalEnv)
    assign("scores_file", scores_file_out, envir = .GlobalEnv)
    assign("unique_infor", unique_infor, envir = .GlobalEnv)
    assign("q_range", q_range, envir = .GlobalEnv)
    assign("q_range_file", q_range_out, envir = .GlobalEnv)
    assign("p_values", p_values, envir = .GlobalEnv)

  }
}
