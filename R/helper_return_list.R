#' Helper function for PreparePlinkFile()
#' @description
#' Return either a list containing the score_file, p_value_file, unique_infor and
#' q_range or the global variables for the same objects and write the files.
#' Not exported
#' @param x return_list
#' @param output path to output folder
#' @param score_file
#' @param p_value_file
#' @param unique_infor
#' @param q_range
#' @return CreateQRange() object
#'
helper_return_list <- function(x,
                               output,
                               score_file,
                               p_value_file,
                               unique_infor,
                               q_range)
  {
  file1 <- paste0(output,"score_file")

  if (x) {
    print("return plink_file")
    result <- list(score_file,
                   p_value_file,
                   unique_infor,
                   q_range)
    assign("plink_files", result, envir = .GlobalEnv)
  } else {
    score_file_out <- paste0(output, "score_file")
    print(paste0("printing score file: ", score_file_out))
    write.table(score_file,
               file = paste0(output,"score_file"),
               row.names = F,
               col.names = F,
               quote=F)
    q_range_out <- paste0(output, "q_range_file")
    print(paste0("printing q_range_file: ", q_range_out))
    write.table(q_range,
               file = paste0(output,"q_range_file"),
               row.names = F,
               col.names = F,
               quote=F)
    assign("score_file", score_file, envir = .GlobalEnv)
    assign("unique_infor", unique_infor, envir = .GlobalEnv)
    assign("q_range", q_range, envir = .GlobalEnv)
    assign("p_value_file", p_value_file, envir = .GlobalEnv)

  }
}
