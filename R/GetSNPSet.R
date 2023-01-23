#' Get the best performance SNPS
#' @description
#' Get the SNP set for estimating the covariance matrix of prior distribution in
#' empirical Bayes method
#' @param snp_set_ind Column name from prs_mat (output from PRSscore()) with best
#' performance against the tuning data set.
#' @param scores The scores dataframe from PreparePlinkFile output.
#' @param clump_info The unique_infor from PreparePlinkFile output.
#'
#' @return data.frame of SNP performance metrics
#' @export
#' @usage GetSNPSe(snp_id, scores, clump_info)

GetSNPSet <- function(snp_ind,
                      scores,
                      clump_info){
  print("Executing GetSNPSet()... ")

  str_temp <- strsplit(snp_ind,"_")
  r2 <- str_temp[[1]][3]
  ws <- str_temp[[1]][5]
  p_ref_cutoff <- as.numeric(str_temp[[1]][[8]])
  p_tar_cutoff <- as.numeric(str_temp[[1]][[11]])

  scores_name <- paste0("clump_r2_",r2,"_ws_",ws)
  idx <- which(colnames(scores)==scores_name)
  ld <- scores[scores[,idx]!=0,"SNP",drop=F]
  ld_infor <- left_join(ld,clump_info,by="SNP")
  snp_ind <- which(ld_infor$P_ref <= p_ref_cutoff|
                     ld_infor$P <= p_tar_cutoff)
  this_snp <- ld_infor[snp_ind,]

  return(this_snp)
}

