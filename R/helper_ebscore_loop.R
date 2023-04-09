#' helper_ebscore_loop
#'
#' Helper for function for PRSscoreEbayes()
#' @param plink2_exec description
#' @param bfile description
#' @param eb_plink_list description
#' @param plink_list description
#' @param pthres Default c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#' @param threads Default 2
#' @param memory Default 8000
#' @param out_prefix description
#' @param params_farm description
#' @param results_dir description

#' @keywords plink1.9 clump
#' @export

helper_ebscore_loop <- function(plink2_exec,
                                bfile,
                                eb_plink_list,
                                plink_list,
                                pthres,
                                threads,
                                memory,
                                results_dir,
                                out_prefix = as.null(),
                                params_farm = as.null()){

  outfile_prefix <- helper_prefix(out_prefix = out_prefix,
                                  ebayes = TRUE)

  scores <- eb_plink_list[[1]]
  p_values_eb <- eb_plink_list[[2]]
  clump_info <- plink_list[[3]]
  #p_values_eb_file <- as.character(unlist(eb_plink_list["p_values_eb_file"]))
  #score_eb_file <- as.character(unlist(eb_plink_list["score_eb_file"]))
  #q_range_file <- as.character(unlist(plink_list["q_range_file"]))
  p_values_eb_file <- paste0(results_dir, "p_values_eb_file")
  score_eb_file <- paste0(results_dir, "score_eb_file")
  q_range_file <- paste0(results_dir, "q_range_file")

  prs_p_other_ <- paste0(results_dir, outfile_prefix, "prs_p_other_")
  assign("prs_p_other_", prs_p_other_, envir = .GlobalEnv)
  p_values_temp <- p_values_eb

  for(k1 in 1:length(pthres)){
    idx <- which(clump_info$P_ref <= pthres[k1])
   # print(paste0("pthres: ", pthres[k1]))
    #print(class(pthres[k1]))
    p_values_temp$P[idx] <- 0
    # print(paste0("Number of variants less than pthres :",
    #              sum(p_values_temp$P == 0)))
    # print(paste0("writing ", p_values_eb_file))
    write.table(p_values_temp,
                file = p_values_eb_file,
                col.names = F,
                row.names = F,
                quote=F)
    score_col_nums <- ncol(scores)
    plink2score(params_farm = params_farm,
                plink2_exec = plink2_exec,
                bfile = bfile,
                q_range_file = q_range_file,
                p_value_file = p_values_eb_file,
                score_file = score_eb_file,
                score_col_nums = score_col_nums,
                results_dir = results_dir,
                pthres_idx = k1,
                threads = threads,
                out = prs_p_other_,
                memory = memory)

  }
  return(prs_p_other_)
}
