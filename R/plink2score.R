#' Generate PRS scores with Plink2
#' @description
#' Performs PRS scoring with plink2 for all combinations of clumpings and p-value
#' threshold combinations
#' @param plink2_exec Plink binary execute file. Default "plink2".
#' @param bfile prefix for target plink binary set (prefix.bed, prefix.bim, prefix.fam)
#' @param q_range_file File name and location of written q_range dataframe produced
#' by PreparePlinkfile(). By default PreparePlinkfile() names this 'q_range_file'
#' @param p_value_file File name and location of written p_values_temp dataframe
#' either produced by helper_score_loop() or manually as indicated in the
#' vignette section 'Generate plink2 PRS'
#' @param score_col_nums number of columns in scores dataframe produced by
#' PreparePlinkfile().
#' @param scores_file File name and location of written scores dataframe produced
#' by PreparePlinkfile(). By default PreparePlinkfile() names this 'scores_file'
#' @param p_values p_values dataframe produced by PreparePlinkfile()
#' @param threads maximum number of concurrent threads
#' @param memory primary workspace memory
#' @param out output folder and file prefix for Plink2. Default name prs_p_other_
#' <pthres column number>.<range name>.sscore
#' @param prs_farm prs_farm object name
#' @keywords plink2 score
#' @usage Plink19Clump(plink2_exec, bfile, q_score_range, score_col_nums,
#' score, threads, memory, out, score_farm)
#' @return plink2 score files
#' @export
#'
plink2score <- function(plink2_exec = "plink2 ",
                        bfile,
                        p_value_file = p_value_file,
                        q_range_file = q_range_file,
                        score_col_nums,
                        scores_file = scores_file,
                        threads = 4,
                        memory = 8000,
                        out,
                        params_farm=as.null())
  {

  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
    mem <- as.character(unlist(score_farm["mem"]))
    threads <- as.character(unlist(score_farm["threads"]))
  }
  system(paste0(plink2_exec, " ",
                "--bfile ", bfile, " ",
                "--q-score-range ", q_range_file, " ", p_value_file, " ",
                "--score-col-nums 3-", score_col_nums, " ",
                "--score ", scores_file, " cols=+scoresums,-scoreavgs ",
                "--threads ", threads, " ",
                "--memory ", memory, " ",
                "--out ", out)
  )
}