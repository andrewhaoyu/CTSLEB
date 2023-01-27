#' Plink --score wrapper
#' @description
#' Wrapper to run plink2a --score routine with EBayes estimates and then combine
#' all the PRSs.
#' @param plink2_exec Plink2 binary execute file. Default "plink2".
#' @param bfile Plink binary files for the target test population.
#' @param eb_plink_list List of plink data frame and files produced from
#' PreparePlinkFileEBayes()
#' @param plink_list List of plink data.frame and files from PreparePlinkFiles()
#' @param pthres P value thresholds. Defaults to
#' c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#' @param threads Plink "--threads" parameter. Defaults to 2
#' @param memory Plink "--memory" parameter. Defaults to 8000 mb
#' @param out_prefix Prefix for exported files. Recommended when plink files
#' are divided by chromosome or chunks, i.e., "chr1" or "chr2", and CTSLEB
#' script is dispatched in parallel.
#' @param params_farm List of plink parameters produced from SetParamsFarm()
#' @param results_dir Folder or directory to export files.
#' @keywords plink2
#' @export
#' @usage PRSscoreEBayes <- function(plink2_exec,bfile,ebayes, eb_plink_list,
#' plink_list, pthres, threads, memory, results_dir, out_prefix, params_farm)
#'
PRSscoreEBayes <- function(plink2_exec = "plink2 ",
                           bfile,
                           eb_plink_list,
                           plink_list,
                           pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                           threads = 2,
                           memory = 8000,
                           results_dir,
                           out_prefix = as.null(),
                           params_farm=as.null()){

  print("executing PRSscoreEbayes()... ")
  if (is.null(params_farm)) {
    #print("no params_farm")
  } else {
    #print("params_farm list will be used")
    plink2_exec <- as.character(unlist(params_farm["plink2_exec"]))
    memory <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }

  this_prs_p_other_ <- helper_ebscore_loop(params_farm = params_farm,
                                           plink2_exec = plink2_exec,
                                           bfile = bfile,
                                           eb_plink_list = eb_plink_list,
                                           plink_list = plink_list,
                                           pthres = pthres,
                                           threads = threads,
                                           memory = memory,
                                           results_dir = results_dir,
                                           out_prefix = out_prefix)
  scores_eb <- eb_plink_list[[1]]
  this_prs_mat <- helper_CombinePRS(scores = scores_eb,
                                    pthres = pthres,
                                    prs_p_other_ = this_prs_p_other_)
  print("prs_mat_eb object created")
  return(this_prs_mat)
}
