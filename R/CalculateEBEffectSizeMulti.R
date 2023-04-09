#' CalculateEBEffectSize
#'
#' Performs the entire Empirical-Bayes estimation of effect sizes
#' workflow as outlined in multi-ancestry section of the vignette
#' @param sum_target The GWAS summary statistics for the target population. The data needs to have following columns at least: CHR, SNP, BP, A1, BETA, SE, P. A1 is the effect allele. BETA is the regression coefficients for linear regression, log-odds ratio for logistic regression. SE is the standard error for BETA.
#' @param sum_ref_list A list that contains the GWAS summary statistics for multiple populations.
#' @param ref_names The names of the other anescestries for GWAS summary statistics
#' @param best_snps_set, Object created by the GetSNPSet() function.
#' @param snp_list, Object created by the RunClump() function.
#' @param bfile Plink binary files for the target test cohort.
#' @param plink2_exec Plink2 binary execute file. Default "plink2".
#' @param out_prefix Prefix for exported files. Recommended when plink files
#' are divided by chromosome or chunks, i.e., "chr1" or "chr2", and CTSLEB
#' script is dispatched in parallel.
#' @param plink_list List of plink data.frame and files from PreparePlinkFiles()
#' @param threads Plink "--threads" parameter. Defaults to 2
#' @param memory Plink "--memory" parameter. Defaults to 8000 mb
#' @param out_prefix Prefix for exported files. Recommended when plink files
#' are divided by chromosome or chunks, i.e., "chr1" or "chr2", and CTSLEB
#' script is dispatched in parallel.
#' @param params_farm List of plink parameters produced from SetParamsFarm()
#' @param results_dir Folder or directory to export files.
#' @keywords EBayes
#' @export
#' @usage CalculateEBEffectSize(bfile,prs_tune,plink_list,memory,threads,
#' out_prefix,results_dir,params_farm)

CalculateEBEffectSizeMulti <- function(ref_names,
                                       sum_ref_list,
                                       bfile,
                                       sum_tar,
                                       best_snps_set,
                                       snp_list,
                                       plink_list,
                                       plink_list_eb,
                                       results_dir,
                                       out_prefix,
                                       threads = 2,
                                       memory = 8000,
                                       params_farm = as.null()){

  this_q_range <- plink_list[[4]]
  clump_info <- plink_list[[3]]  # (i.e., unique_infor)

  multi_sum_com <- AlignSumMulti(sum_tar = sum_AFR,
                                 sum_ref_list = sum_other_list,
                                 ref_names = ref_names)
  assign("multi_sum_com", multi_sum_com, envir = .GlobalEnv)

  multi_unique_infor_post <- EBpostMulti(x = clump_info,
                                         y = best_snps_set,
                                         sum_com = multi_sum_com,
                                         ref_names = ref_names)
  assign("multi_unique_infor_post", multi_unique_infor_post, envir = .GlobalEnv)

  multi_eb_post_col_names <- c("BETA_EB_target",paste0("BETA_EB_",ref_names[1]))
  multi_post_beta_mat <- multi_unique_infor_post %>%
    select(all_of(multi_eb_post_col_names))

  multi_plink_list_eb <- PreparePlinkFileEBayes(snp_list = snp_list,
                                               clump_info = clump_info,
                                               post_clump_info = multi_unique_infor_post,
                                               post_beta = multi_post_beta_mat,
                                               results_dir = results_dir)
  assign("multi_plink_list_eb", multi_plink_list_eb, envir = .GlobalEnv)

  multi_score_eb <- multi_plink_list_eb[['scores_eb']]
  assign("multi_score_eb", multi_score_eb, envir = .GlobalEnv)

  score_eb_file <- as.character(unlist(multi_plink_list_eb["score_eb_file"]))
  write.table(multi_score_eb,
              file = score_eb_file,
              row.names = F,
              col.names = F,
              quote=F)

  write.table(this_q_range,
              file = paste0(results_dir,"q_range_file"),
              row.names = F,
              col.names = F,
              quote=F)
  #TODO:  check id clump_info (i.e., unique_infor) can be passed instead
  # of plink_list
  multi_prs_mat_eb <- PRSscoreEBayes(bfile = bfile,
                                    eb_plink_list = multi_plink_list_eb,
                                    plink_list = plink_list,
                                    results_dir = results_dir,
                                    out_prefix = out_prefix,
                                    threads = threads,
                                    memory = memory,
                                    params_farm = PRS_farm)

  colnames(multi_prs_mat_eb)[1] <- "FID"
  return(multi_prs_mat_eb)
}

