#' Prepare the data.frames and file names for the PLINK2 --score routine using
#' EBayes coefficients
#' @description Create the PRS coefficients matrix for all ancestries.
#' @param snp_list description
#' @param results_dir description
#' @param clump_info description
#' @param post_clump_info description
#' @param post_beta description
#' @keywords plink1.9 clump
#' @export
#' @usage PreparePlinkFileEBayes(snp_list, clump_info, post_clump_info,
#' post_beta, results_dir)
#' @return List of data.frames and file locations for Plink2; scores_eb [[1]],
#' p_values_eb [[2]], score_eb_file [[3]], p_values_eb_file. The scores_eb data.frame
#' first column contains the unique SNPs after LD clumping. The second column is
#' the effect allele. The third through last columns contain the EBayes posterior
#' mean for the reference and target populations for SNPs after specific LD clumping
#' combinations. The coefficients are set to 0 if a SNP doesn't exist for the
#' clumping combination. The p_value_eb data.frame first column is the same as
#' scores_eb. The second column is the p-values for SNPs from the target population GWAS.

PreparePlinkFileEBayes <- function(snp_list,
                                   clump_info,
                                   post_clump_info,
                                   post_beta,
                                   results_dir){

  print("Executing PreparePlinkFileEBayes() ... ")
  unique_id <- post_clump_info$SNP
  names(unique_id) <- "SNP"

  #number of ancestry
  n_ans <- ncol(post_beta)

  n_col = length(snp_list)
  n_row <- nrow(post_clump_info)
  post_beta <- as.matrix(post_beta)
  post_beta[is.na(post_beta)] <- 0

  beta_mat <- matrix(rep(post_beta,n_col),nrow =n_row,ncol =n_col*n_ans)
  names <- rep("c",n_col*n_ans)
  temp <- 0
  for(ldx in 1:n_col){
    ld <- snp_list[[ldx]]
    names(ld) = "SNP"
    idx <- which(clump_info$SNP%in%ld$SNP==F)
    beta_mat[idx,(1:n_ans)+temp] = 0
    names[(1:n_ans)+temp] <- paste0(names(snp_list[[ldx]]),
                                    "_",
                                    colnames(post_beta))
    temp <- temp + n_ans
  }
  colnames(beta_mat) <- names
  this_scores <- data.frame(SNP = unique_id,
                            A1 = post_clump_info$A1,
                            beta_mat)
  temp_dir <- paste0(results_dir,"temp/")
  score_eb_file <- paste0(temp_dir,"score_eb_file")
  this_p_values <- data.frame(SNP = unique_id,
                              P = post_clump_info$P)
  p_values_eb_file <-paste0(temp_dir,"p_values_eb_file")
  names <- c("scores_eb",
             "p_values_eb",
             "score_eb_file",
             "p_values_eb_file")
  values <- list(this_scores,
                 this_p_values,
                 score_eb_file,
                 p_values_eb_file)

  print("PreparePlinkFileEBayes() complete ...")
  result <- setNames(values, names)

  return(result)
}
