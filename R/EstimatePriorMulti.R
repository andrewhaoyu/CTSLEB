#' Estimate the covariance matrix for the prior distribution used in EB step
#'
#' @param  SNP_set the SNP set from CT step for estimating the covariance matrix for the prior distribution
#'
#' @return the covariance matrix estimate
#'
#' @examples
EstimatePriorMulti <- function(SNP_set,other_ans_names,
                               sum_com){
  SNP_set_select = SNP_set %>%
    select(SNP)
  #align the summary statistics with the SNP_set from CT
  SNP_set_align = left_join(SNP_set_select,sum_com,
                            by="SNP")
  n_ans = length(other_ans_names)

  #create column names to select the BETA_ and SE_ columns
  col_names_beta = c("BETA",paste0("BETA_",other_ans_names))
  col_names_se = c("SE",paste0("SE_",other_ans_names))
  beta_mat = SNP_set_align %>%
    select(all_of(col_names_beta))
  se_mat = SNP_set_align %>%
    select(all_of(col_names_se))
  #\hat_{u}_kl|u_kl ~ N(u_kl,1/N_l)
  #where u_kl is the underlying effect size for the kth SNP of lth population
  #N_l is the sample size
  #since the Bayesian algorithm is applied on the standardized effect-size scale
  #it's equivalent to applying the Bayes rule on z-statistics scale
  #the advantage of z-statistics scale is the covariance matrix is identity
  #it can make the computation faster.
  z_mat = beta_mat/se_mat
  z_mat <-na.omit(z_mat)
  p = ncol(z_mat)
  prior_mat <- cov(z_mat)-diag(p)
  colnames(prior_mat) = c("Z_tar",paste0("Z_",other_ans_names))
  return(prior_mat)
}

