#' Create the q range file used for plink2
#'
#' @param pthres the p value thresholds for selecting SNPs
#'
#' @return the q range file
#' @export
#'
CreateQRange <- function(pthres){
  n_pthres = length(pthres)
  q_range = data.frame(
    filename = rep("p_tar",n_pthres),
    small_P = rep(0,n_pthres),
    max_P = rep(1,n_pthres))

  temp = 1

  for(k2 in 1:length(pthres)){

    q_range[temp,1] = paste0("p_tar_",k2)
    q_range[temp,3] = pthres[k2]
    temp = temp+1
  }
  q_range = q_range[1:(temp-1),]
  return(q_range)
}
