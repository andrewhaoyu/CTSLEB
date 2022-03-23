#' Estimate Bayesian posterior meanusing the empirical Bayes algorithm
#'
#' @param unique_infor the unique_infor from PreparePlinkFile output, which contains information for all SNPs after clumping step
#' @param SNP_set the SNP set from CT step for estimating the covariance matrix for the prior distribution

#' @export
#' @return Bayesian postier mean for all SNPs after-clumping
EBpost <- function(unique_infor,SNP_set){
  EBprior = EstimatePrior(SNP_set)
  BETA_tar = as.numeric(unique_infor$BETA)
  SE_tar = as.numeric(unique_infor$SE)
  BETA_other = as.numeric(unique_infor$BETA_other)
  SE_other = as.numeric(unique_infor$SE_other)
  #\hat_{u}_kl|u_kl ~ N(u_kl,1/N_l)
  #where u_kl is the underlying effect size for the kth SNP of lth population
  #N_l is the sample size
  #since the Bayesian algorithm is applied on the standardized effect-size scale
  #it's equivalent to applying the Bayes rule on z-statistics scale
  #the advantage of z-statistics scale is the covariance matrix is identity
  #it can make the computation faster.
  prior_sigma = EBprior

  z_tar = BETA_tar/SE_tar
  z_other = BETA_other/SE_other
  z_mat <-as.matrix(cbind(z_tar,z_other))

  post_sigma = solve(solve(prior_sigma)+diag(2))
  z_mat_post = z_mat

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
  beta_mat_post = z_mat_post
  beta_mat_post[,1] =z_mat_post[,1]*SE_tar
  beta_mat_post[,2] =z_mat_post[,2]*SE_other
  colnames(beta_mat_post) = c("BETA_EB_target","BETA_EB_other")
  unique_infor_EB =cbind(unique_infor,beta_mat_post) %>%
    select(SNP,A1,BETA_EB_target,P,BETA_EB_other,P_other)
  return(unique_infor_EB)
}
