#' Helper function for dimCT()
#' @description
#' Inputs list object from PreparePlinkFile() and creates global variables
#' and write.tables.
#' q_range from Prepare
#' Not exported
#' @param plink_list description
#' @param results_dir description
#' @return list
#'
helper_PreparePlinkFile <- function(plink_list,
                                    results_dir)

  {

  temp.dir <- paste0(results_dir,"temp/")
  scores <- plink_list[[1]]
  score_file <- paste0(temp.dir,"score_file")
  write.table(scores,
              file = score_file ,
              row.names = F,
              col.names = F,
              quote=F)

  p_values <- plink_list[[2]]
  p_value_file <- paste0(temp.dir,"p_value_file")

  unique_infor <- plink_list[[3]]

  q_range <- plink_list[[4]]
  q_range_file <- paste0(temp.dir,"q_range_file")

  write.table(q_range,
              file = q_range_file,
              row.names = F,
              col.names = F,
              quote=F)
  assign("q_range", q_range, envir = .GlobalEnv)
  assign("scores", scores, envir = .GlobalEnv)
  assign("p_values", p_values, envir = .GlobalEnv)
  assign("unique_infor", unique_infor, envir = .GlobalEnv)
  assign("score_file", score_file, envir = .GlobalEnv)
  assign("p_value_file", p_value_file, envir = .GlobalEnv)
  assign("q_range_file", q_range_file, envir = .GlobalEnv)

  names <-c("score_file","p_value_file","q_range_file")
  values <- list(score_file, p_value_file, q_range_file)
  list <- setNames(values, names)

  return(list)
}
