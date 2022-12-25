#' Helper function for dimCT()
#' @description
#' Inputs list object from PreparePlinkFile() and creates global variables
#' and write.tables.
#' q_range from Prepare
#' Not exported
#' @param plink_list
#' @param results_dir
#' @return
#'
helper_PreparePlinkFile <- function(plink_list,
                                    results_dir)

  {
  this.plink_list <- plink_list
  this.results_dir <- results_dir

  temp.dir <- paste0(results_dir,"temp/")
  scores <- this.plink_list[[1]]
  scores_file <- paste0(temp.dir,"score_file")
  write.table(scores,
              file = scores_file ,
              row.names = F,
              col.names = F,
              quote=F)

  p_values <- this.plink_list[[2]]
  unique_infor <- this.plink_list[[3]]

  q_range <- plink_list[[4]]
  head(q_range)
  q_range_file <- paste0(temp.dir,"q_range_file")
  write.table(q_range,
              file = q_range_file,
              row.names = F,
              col.names = F,
              quote=F)
    assign("scores", scores, envir = .GlobalEnv)
    assign("scores_file", scores_file, envir = .GlobalEnv)
    assign("unique_infor", unique_infor, envir = .GlobalEnv)
    assign("q_range", q_range, envir = .GlobalEnv)
    assign("q_range_file", q_range_file, envir = .GlobalEnv)
    assign("p_values", p_values, envir = .GlobalEnv)
}
