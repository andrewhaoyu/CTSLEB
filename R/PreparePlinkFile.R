#' Prepare the files for PLINK2 to calculate PRSs
#'
#' @param snp_list the snp_list result from the two-dimensional clumping
#' @param sum_com the sum_com result from the AlignSum function.
#'
#' @return the score_file, p_value_file and p_other_file for the two-dimensional thresholding
#' @export
#'
#' @examples
PreparePlinkFile = function(snp_list,
                            sum_com){
  #create unique SNP list by combind LD clumping results under different parameters
  unique_id = unique(rbindlist(snp_list,use.name =FALSE))
  names(unique_id) = "SNP"
  #align the regression coefficients for these SNPs from the sum stat
  unique_infor = left_join(unique_id,sum_com,by="SNP")

  #create a coefficient matrix
  #the first column contains the unique SNPs after clumping results under all combinations of r2-cutoff and window_size
  #the second column is the effect allele
  #the third to the last columns contains the regression coefficients of the target population for SNPs after LD-clumping under a specific combination of r2-cutoff and base_window_size
  #the coefficients is put as 0 if a SNP doesn't exist in the clumping results under a specific combination of r2-cutoff and base_window_size
  n_col = length(snp_list)
  n_row = nrow(unique_infor)
  beta_mat = matrix(unique_infor$BETA,nrow =n_row,ncol =n_col)
  names = rep("c",n_col)
  temp = 1
  for(ldx in 1:n_col){
    LD = snp_list[[ldx]]
    names(LD) = "SNP"
    idx <- which(unique_infor$SNP%in%LD$SNP==F)
    beta_mat[idx,ldx] = 0
    names[ldx] = names(snp_list[[ldx]])
  }
  colnames(beta_mat) = names
  score_file = data.frame(SNP = unique_id,A1 = unique_infor$A1,beta_mat)
  p_value_file = data.frame(SNP = unique_id,P = unique_infor$P)
  result = list(score_file,
                p_value_file,
                unique_infor)
  return(result)
}
