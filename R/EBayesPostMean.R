#' EBayesPostMean
#'
#' This function performs the entire Empirical-Bayes estimation of effect sizes
#' workflow as outlined in step2 of the vignette
#' @param clump_info description
#' @param snp_set description
#' @usage EBayesPostMean(clump_info, snp_set)
#' @export

EBayesPostMean <- function(clump_info,snp_set){
  print("executing EBayesPostMean()... ")
  prior_sigma <- EstimatePrior(snp_set)
  beta_tar <- as.numeric(clump_info$BETA)
  se_tar <- as.numeric(clump_info$SE)
  beta_ref <- as.numeric(clump_info$BETA_ref)
  se_ref <- as.numeric(clump_info$SE_ref)
  z_tar <- beta_tar/se_tar
  z_ref <- beta_ref/se_ref
  z_mat <-as.matrix(cbind(z_tar,z_ref))

  post_sigma <- solve(solve(prior_sigma)+diag(2))
  z_mat_post <- z_mat

  p <- ncol(z_mat)

  for(k in 1:nrow(z_mat)){
    if(k%%10000==0){print(paste0(k," SNPs completed"))}
    z_temp = z_mat[k,]

    #find out nonmissing component

    idx <- which(!is.na(z_temp))
    if(length(idx)<p){
      z_temp <- z_temp[idx]

      post_sigma_temp = post_sigma[idx,idx,drop=F]
      z_post = post_sigma_temp%*%z_temp
    }else{
      z_post =post_sigma%*%z_temp
    }

    z_mat_post[k,idx] = z_post
  }
  beta_mat_post <- z_mat_post
  beta_mat_post[,1]  <- z_mat_post[,1]*se_tar
  beta_mat_post[,2] <- z_mat_post[,2]*se_ref
  colnames(beta_mat_post) <- c("BETA_EB_target","BETA_EB_ref")
  clump_info_ebayes <- cbind(clump_info,beta_mat_post) %>%
    select(SNP,A1,BETA_EB_target,P,BETA_EB_ref,P_ref)
  print("EBpost() complete... ")

  return(clump_info_ebayes)
}
