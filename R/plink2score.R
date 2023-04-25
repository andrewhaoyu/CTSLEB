#' Generate PRS scores with Plink2
#' @description
#' Performs PRS scoring with plink2 for all combinations of clumpings and p-value
#' threshold combinations
#' @param plink2_exec Plink binary execute file. Default "plink2".
#' @param bfile prefix for target plink binary set (prefix.bed, prefix.bim, prefix.fam)
#' @param q_range_file File name and location of written q_range dataframe produced
#' by PreparePlinkfile(). By default PreparePlinkfile() names this 'q_range_file'
#' @param p_value_file description
#' @param scores_file File name and location of written scores dataframe produced
#' by PreparePlinkfile(). By default PreparePlinkfile() names this 'scores_file'
#' @param score_col_nums number of columns in scores dataframe produced by
#' PreparePlinkfile().
#' @param results_dir output folder forvPlink2. Results will appear as results_dir/
#' temp/prs_p_other_<pthres index>.<range name>.sscore
#' @param pthres_idx pthres vector index number
#' @param threads maximum number of concurrent threads
#' @param memory primary workspace memory
#' @param out description
#' @param params_farm prs_farm object name
#' @keywords plink2 score
#' @usage Plink19Clump(plink2_exec, bfile, q_score_range, score_col_nums,
#' score, threads, memory, out, score_farm)
#' @return plink2 score files
#' @export
#'
plink2score <- function(plink2_exec = "plink2 ",
                        bfile,
                        q_range_file,
                        p_value_file,
                        score_file,
                        score_col_nums,
                        results_dir,
                        pthres_idx,
                        threads = 2,
                        memory = 8000,
                        out,
                        params_farm=as.null())
  {

  if (is.null(params_farm)) {
    #print("no params_farm")
  } else {
    #print("params_farm list will be used")
    plink2_exec <- as.character(unlist(params_farm["plink2_exec"]))
    mem <- as.character(unlist(params_farm["mem"]))
    threads <- as.character(unlist(params_farm["threads"]))
  }

  prs_out <- paste0(prs_p_other_,pthres_idx)
  If(score_col_nums==3){

    system(paste0(plink2_exec, " ",
                  "--bfile ", bfile, " ",
                  "--q-score-range ",  q_range_file, " ", p_value_file, " ",
                  "--score-col-nums ",  score_col_nums, " ",
                  "--score ", score_file, " cols=+scoresums,-scoreavgs ",
                  "--threads ", threads, " ",
                  "--memory ", memory, " ",
                  "--out ", prs_out))

  }else{
    system(paste0(plink2_exec, " ",
                  "--bfile ", bfile, " ",
                  "--q-score-range ", q_range_file, " ", p_value_file, " ",
                  "--score-col-nums 3-", score_col_nums, " ",
                  "--score ", score_file, " cols=+scoresums,-scoreavgs ",
                  "--threads ", threads, " ",
                  "--memory ", memory, " ",
                  "--out ", prs_out))
  }

  )
}
