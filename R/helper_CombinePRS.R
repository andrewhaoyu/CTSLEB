#' Helper function for PRSscore()
#' @description
#' Combines PRS
#' Not exported
#' @param scores  description
#' @param pthres description
#' @param prs_p_other_  description
#' @return 'prs_mat' global variable. Matrix of combined PRSs with clumping and
#' pthres as column name
#'
helper_CombinePRS <- function(scores,
                              pthres,
                              prs_p_other_)
  {

  print("combining PRS results ...")
  prs_list <- list()
  temp <- 1
  names <- colnames(scores[3:ncol(scores)])
  for(k1 in 1:length(pthres)){
    for(k2 in 1:length(pthres)){
      prs_file <- paste0(prs_p_other_,k1,".p_tar_",k2,".sscore")
      #print(paste0("reading :", prs_file))
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
  print("prs_mat complete... ")
  return(prs_mat)
}
