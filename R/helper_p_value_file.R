#' Helper function for PreparePlinkFile()
#' @description
#' Find reference SNPs with p-value less than pthres and set the p-value to
#' 0 to guarantee their selection.
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
    idx <- which(unique_infor$P_ref<=pthres[k1])
    p_value_file_temp$P[idx] <- 0
    write.table(p_value_file_temp,
                file = paste0(output,"p_value_file"),
                col.names = F,
                row.names = F,
                quote=F)
  }

}
