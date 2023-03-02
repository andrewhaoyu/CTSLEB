#' Estimate Bayesian posterior meanusing the empirical Bayes algorithm
#'
#' @param x unique_infor the unique_infor from PreparePlinkFile output, which contains information for all SNPs after clumping step
#' @param y SNP_set the SNP set from CT step for estimating the covariance matrix for the prior distribution
#' @param sum_com summary statistics for all ancestries after allele alignment
#' @param ref_names the other ancestries' name used in the analyses
#' @export
#' @return Bayesian postier mean for all SNPs after-clumping

EBpostMulti <- function(x,
                        y,
                        sum_com,
                        ref_names){
  print("executing EBpostMulti()... ")
  prior_sigma <- EstimatePriorMulti(snp_set = y,
                                    ref_names = ref_names,
                                    sum_com = sum_com)
  SNP_set_select <- x %>% select(SNP)

  #align the summary statistics with the SNP_set from CT

  SNP_set_align <- left_join(SNP_set_select,sum_com,
                             by="SNP")


  #create column names to select the BETA_ and SE_ columns

  col_names_beta <- c("BETA",paste0("BETA_",ref_names))
  col_names_se <- c("SE",paste0("SE_",ref_names))
  beta_mat <- SNP_set_align %>% select(all_of(col_names_beta))
  se_mat <- SNP_set_align %>% select(all_of(col_names_se))

  z_mat <- as.matrix(beta_mat/se_mat)
  z_mat_post <- as.matrix(z_mat)
  col_names_beta <- c("Z",paste0("Z_",ref_names))
  p <- ncol(z_mat)

  post_sigma <- solve(solve(prior_sigma)+diag(p))


  for(k in 1:nrow(z_mat)){
    if(k%%10000==0){print(paste0(k," SNPs completed"))}
    z_temp <- z_mat[k,]

    #find out nonmissing component

    idx <- which(!is.na(z_temp))
    if(length(idx)<p){
      z_temp <- z_temp[idx]

      post_sigma_temp <- post_sigma[idx,idx,drop=F]
      z_post <- post_sigma_temp%*%z_temp
    }else{
      z_post <- post_sigma%*%z_temp
    }

    z_mat_post[k,idx] <- z_post
  }
  beta_mat_post <- z_mat_post*se_mat
  colnames(beta_mat_post) <- c("BETA_EB_target",paste0("BETA_EB_",ref_names))
  eb_beta_names <- colnames(beta_mat_post)
  unique_infor_EB <- cbind(unique_infor,beta_mat_post) %>%
    select(SNP,A1,all_of(eb_beta_names),P,P_ref)
  return(unique_infor_EB)
}
