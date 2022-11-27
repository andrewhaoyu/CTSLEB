#' Helper function for PreparePlinkFile()
#' @description
#' Either a list containing the score_file, p_value_file, unique_infor and
#' q_range or the global variables for the same objects and writes the files.
#' Not exported
#' @param x p_value_file
#' @return updated p_value_file object
#'
helper_p_value_file <- function(x,pthres,output)
  {
  #create a temporary p_value_file
  p_value_file_temp <- x
  for(k1 in 1:length(pthres)){
    #keep all the SNPs with P_EUR less than pthres[k1] in the analyses
    idx <- which(unique_infor$P_other<=pthres[k1])
    p_value_file_temp$P[idx] <- 0
    write.table(p_value_file_temp,
                file = paste0(output,"p_value_file"),
                col.names = F,
                row.names = F,
                quote=F)
    n_col <- ncol(score_file)
  }
}
