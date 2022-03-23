#' Estimate the covariance matrix for the prior distribution used in EB step
#'
#' @param  SNP_set the SNP set from CT step for estimating the covariance matrix for the prior distribution
#'
#' @return the covariance matrix estimate
#'
#' @examples
EstimatePrior <- function(SNP_set){
  BETA_tar = as.numeric(SNP_set$BETA)
  SE_tar = as.numeric(SNP_set$SE)
  BETA_other = as.numeric(SNP_set$BETA_other)
  SE_other = as.numeric(SNP_set$SE_other)

  #\hat_{u}_kl|u_kl ~ N(u_kl,1/N_l)
  #where u_kl is the underlying effect size for the kth SNP of lth population
  #N_l is the sample size
  #since the Bayesian algorithm is applied on the standardized effect-size scale
  #it's equivalent to applying the Bayes rule on z-statistics scale
  #the advantage of z-statistics scale is the covariance matrix is identity
  #it can make the computation faster.
  z_tar = BETA_tar/SE_tar
  z_other = BETA_other/SE_other

  z_mat <-na.omit(cbind(z_tar,z_other))


  prior.mat <- cov(z_mat)-diag(2)
  return(prior.mat)
}
