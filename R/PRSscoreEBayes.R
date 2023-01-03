#' PRSscoreEbayes
#'
#' Performs plink2a --score with EBayes estimates
#' @param plink2_exec description
#' @param bfile description
#' @param eb_plink_list description
#' @param plink_list description
#' @param pthres Default c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#' @param threads Default 4
#' @param memory Default 8000
#' @param out_prefix description
#' @param params_farm description
#' @param results_dir description

#' @keywords plink1.9 clump
#' @export

PRSscoreEBayes <- function(plink2_exec = "plink2 ",
                           bfile,
                           ebayes = TRUE,
                           eb_plink_list,
                           plink_list,
                           pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                           threads = 4,
                           memory = 8000,
                           results_dir,
                           out_prefix = as.null(),
                           params_farm=as.null()){

  print("executing PRSscore()... ")
  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
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
