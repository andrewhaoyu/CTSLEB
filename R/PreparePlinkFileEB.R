#' Prepare the files for PLINK2 to calculate PRSs
#'
#' @param snp_list the snp_list result from the two-dimensional clumping
#' @param sum_com the sum_com result from the AlignSum function.
#'
#' @return the score_file, p_value_file and p_other_file for the two-dimensional thresholding
#' @export
#'
#' @examples
PreparePlinkFileEB = function(snp_list,
                            unique_infor_post,
                            post_beta_mat){
  #create unique SNP list by combind LD clumping results under different parameters
  unique_id = unique_infor_post$SNP
  names(unique_id) = "SNP"

  #create a coefficient matrix
  #the first column contains the unique SNPs after clumping results under all combinations of r2-cutoff and window_size
  #the second column is the effect allele
  #the third to the last columns contains the EB coefficients of both the target and EUR population for SNPs after LD-clumping under a specific combination of r2-cutoff and base_window_size
  #the coefficients is put as 0 if a SNP doesn't exist in the clumping results under a specific combination of r2-cutoff and base_window_size
  n_col = length(snp_list)
  n_row = nrow(unique_infor_post)

  #number of ancestry
  n_ans = ncol(post_beta_mat)
  post_beta_mat = as.matrix(post_beta_mat)
  post_beta_mat[is.na(post_beta_mat)] = 0

  beta_mat = matrix(rep(post_beta_mat,n_col),nrow =n_row,ncol =n_col*n_ans)
  names = rep("c",n_col*n_ans)
  temp = 0
  for(ldx in 1:n_col){
    LD = snp_list[[ldx]]
    names(LD) = "SNP"
    idx <- which(unique_infor$SNP%in%LD$SNP==F)
    beta_mat[idx,(1:n_ans)+temp] = 0
    names[(1:n_ans)+temp] = paste0(names(snp_list[[ldx]]),"_",colnames(post_beta_mat))
    temp = temp + n_ans
  }
  colnames(beta_mat) = names
  score_file = data.frame(SNP = unique_id,A1 = unique_infor_post$A1,beta_mat)
  p_value_file = data.frame(SNP = unique_id,P = unique_infor_post$P)
  result = list(score_file,
                p_value_file)
  return(result)
}
