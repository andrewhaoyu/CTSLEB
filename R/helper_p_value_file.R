#' Helper function for PreparePlinkFile()
#' @description
#' Find reference SNPs with p-value less than pthres and set the p-value to
#' 0 and guarantee their selection.
#' Not exported
#' @param x p_value_file_temp
#' @return updated p_value_file object
#'
helper_p_value_file <- function(x = p_value_file,
                                pthres = pthres,
                                unique_infor = unique_infor,
                                score_file = score_file,
                                output = "./")
  {
  #create a temporary p_value_file
  p_value_file <- x
  for(k1 in 1:length(pthres)){
    #keep all the SNPs with P_EUR less than pthres[k1] in the analyses
    idx <- which(unique_infor$P_ref<=pthres[k1])
    p_value_file$P[idx] <- 0
    write.table(p_value_file,
                file = paste0(output,"p_value_file"),
                col.names = F,
                row.names = F,
                quote=F)
    n_col <- ncol(score_file)
  }
  return(p_value_file)

}
