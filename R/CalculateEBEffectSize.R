#' CalculateEBEffectSize
#'
#' Performs the entire Empirical-Bayes estimation of effect sizes
#' workflow as outlined in step2 of the vignette
#' @param bfile Plink binary files for the target test population.
#' @param plink2_exec Plink2 binary execute file. Default "plink2".
#' @param snp_ind Column name from prs_mat (output from PRSscore()) with best
#' performance against the tuning data set.
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

CalculateEBEffectSize <- function(bfile,
                                  snp_ind,
                                  plink_list,
                                  memory = 8000,
                                  threads = 2,
                                  out_prefix,
                                  results_dir,
                                  params_farm = as.null()){
  print("Executing CalculateEBEffectSize() ... ")
  scores <- plink_list[[1]]
  clump_info <- plink_list[[3]]

  best_snp_set <- GetSNPSet(snp_ind = snp_ind,
                            scores = scores,
                            clump_info = clump_info)

  clump_info_post <- EBayesPostMean(x = clump_info,
                                    y = best_snp_set)
  post_coef_mat <- cbind(clump_info_post$BETA_EB_target,
                         clump_info_post$BETA_EB_ref)
  colnames(post_coef_mat) <- c("EB_target","EB_ref")

  plinklist_eb <- PreparePlinkFileEBayes(snp_list = snp_list,
                                          clump_info = clump_info,
                                          post_clump_info = clump_info_post,
                                          post_beta = post_coef_mat,
                                          results_dir = results_dir)

  assign("best_snps_set", best_snp_set, envir = .GlobalEnv)
  assign("unique_infor_post", clump_info_post, envir = .GlobalEnv)
  assign("plink_list_eb", plinklist_eb, envir = .GlobalEnv)

  scores_eb <- plinklist_eb[[1]]
  score_eb_file <- as.character(unlist(plinklist_eb["score_eb_file"]))
  write.table(scores_eb,
              file = score_eb_file,
              row.names = F,
              col.names = F,
              quote=F)
  p_values_eb <- plinklist_eb[[2]]

  ebayes_prs <- PRSscoreEBayes(bfile = bfile,
                               eb_plink_list = plinklist_eb,
                               plink_list = plink_list,
                               results_dir = results_dir,
                               out_prefix = out_prefix,
                               params_farm = params_farm)

  assign("scores_eb", scores_eb, envir = .GlobalEnv)
  assign("score_eb_file", score_eb_file, envir = .GlobalEnv)
  assign("p_values_eb", p_values_eb, envir = .GlobalEnv)

  return(ebayes_prs)

}

