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

# GetSNPSet = function(snp_set_ind,
#                      score_file,
#                     unique_infor){
#   str_temp = strsplit(snp_set_ind,"_")
#   r2 = str_temp[[1]][3]
#   ws = str_temp[[1]][5]
#   p_other_cutoff = as.numeric(str_temp[[1]][[8]])
#   p_tar_cutoff = as.numeric(str_temp[[1]][[11]])
#
#   score_file_name = paste0("clump_r2_",r2,"_ws_",ws)
#   idx = which(colnames(score_file)==score_file_name)
#
#   #take the post-clupmping SNPs with particular r2-cutoff and window size by removing SNPs with coefficients 0
#   LD = score_file[score_file[,idx]!=0,"SNP",drop=F]
#   #take the p_value cutoff
#   LD_infor = left_join(LD,unique_infor,by="SNP")
#   #take the SNP with p_eur < p_other_cutoff | p_tar < p_tar_cutoff
#   snp_ind = which(LD_infor$P_other<=p_other_cutoff|
#                     LD_infor$P<=p_tar_cutoff)
#   SNP = LD_infor[snp_ind,]
#   return(SNP)
#}
