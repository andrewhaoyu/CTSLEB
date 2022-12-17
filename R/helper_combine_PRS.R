#' Helper function for PRSscore()
#' @description
#' Combines PRS
#' Not exported
#' @param scores
#' @param pthres
#' @param prs_file_prefix
#' @return 'prs_mat' global variable. Matrix of combined PRSs with clumping and
#' pthres as column name
#'
helper_combine_PRS <- function(scores = scores,
                               pthres = pthres,
                               prs_file_prefix)
  {
  prs_list <- list()
  temp <- 1
  names <- colnames(scores[3:ncol(scores)])
  for(k1 in 1:length(pthres)){
    for(k2 in 1:length(pthres)){
      prs_file <- paste0(prs_file_prefix,"_",k1,".p_tar_",k2,".sscore")
      prs_temp <- fread(prs_file)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      colnames(prs_list[[temp]]) <- paste0(names,"_",
                                           "p_other_",
                                           pthres[k1],
                                           "_p_tar_",
                                           pthres[k2])
      temp <- temp + 1
    }
  }
  prs_mat <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
  return(prs_mat)
}
