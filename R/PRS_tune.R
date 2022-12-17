#' Select best set of snps for PRS
#' @param x Plink binary execute file. Default "plink2".
#' @param r2_vec prefix for target plink binary set (prefix.bed, prefix.bim, prefix.fam)
#' @param pthres File name and location of written q_range dataframe produced
#' by PreparePlinkfile(). By default PreparePlinkfile() names this 'q_range_file'
#' @param wc_base_vec scores dataframe produced by PreparePlinkfile().
#' @param pheno_df p_values dataframe produced by PreparePlinkfile()
#' @param tune_number Number of samples to select from x for tuning
#' @keywords PRS, tuning
#' @usage Plink19Clump(plink2_exec, bfile, q_score_range, score_col_nums,
#' score, threads, memory, out, score_farm)
#' @export
#' @examples


PRS_tune <- function(x,
                     r2_vec,
                     pthres,
                     wc_base_vec,
                     pheno_df,
                     tune_samples,
                     out = "prs_p_other",
                     params_farm=as.null()){


  prs_mat <- x
  y_tune <- pheno_df
  prs_tune <- prs_mat[1:tune_samples,]
  n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
  prs_r2_vec_test <- rep(0,n.total.prs)
  for(p_ind in 1:n.total.prs){
    #the first two columns of prs_tun are family id and individual id
    #prs starts from the third column
    model <- lm(y_tune$V1~prs_tune[,(2+p_ind)])
    prs_r2_vec_test[p_ind] <- summary(model)$r.square
  }
  max_ind <- which.max(prs_r2_vec_test)
  #+2 is due to the first two columns are family id and individual id
  print(colnames(prs_tune)[max_ind+2])
}
