#' Estimate the covariance matrix for the prior distribution used in EB step
#'
#' @param  snp_set the SNP set from CT step for estimating the covariance matrix for the prior distribution
#' @param  ref_names description
#' @param  sum_com AlignSumMulti() object
#' @return the covariance matrix estimate
#' @export

EstimatePriorMulti <- function(snp_set,
                               ref_names,
                               sum_com){
  print("executing EstimatePriorMulti()... ")
  SNP_set_select <- snp_set %>% select(SNP)

  #align the summary statistics with the SNP_set from CT

  SNP_set_align <- left_join(SNP_set_select,
                             sum_com,
                             by="SNP")
  n_ans <- length(ref_names)

  #create column names to select the BETA_ and SE_ columns

  col_names_beta <- c("BETA",paste0("BETA_",ref_names))
  col_names_se <- c("SE",paste0("SE_",ref_names))
  beta_mat <- SNP_set_align %>% select(all_of(col_names_beta))
  se_mat <- SNP_set_align %>% select(all_of(col_names_se))
  #\hat_{u}_kl|u_kl ~ N(u_kl,1/N_l)
  #where u_kl is the underlying effect size for the kth SNP of lth population
  #N_l is the sample size
  #since the Bayesian algorithm is applied on the standardized effect-size scale
  #it's equivalent to applying the Bayes rule on z-statistics scale
  #the advantage of z-statistics scale is the covariance matrix is identity
  #it can make the computation faster.
  z_mat <- beta_mat/se_mat
  z_mat <-na.omit(z_mat)
  p <- ncol(z_mat)
  prior_mat <- cov(z_mat)-diag(p)
  colnames(prior_mat) <- c("Z_tar",paste0("Z_",ref_names))
  print("EstimatePriorMulti() complete... ")
  return(prior_mat)
}
