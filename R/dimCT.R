#' dimCT
#'
#' This function performs the entire two-dimensional clumping and thresholding
#' workflow (dimCT) as outlined in step1 of the vignette
#' @param plink19_exec Plink binary execute file. Default "plink".
#' @param plink2_exec Plink2 binary execute file. Default "plink2".
#' @param results_dir Folder or directory to export files
#' @param sum_ref The GWAS summary statistics for the reference population.
#' @param sum_target The GWAS summary statistics for the target population.
#' @param ref_plink Plink binary files for the reference population.
#' @param target_plink Plink binary files for the target population.
#' @param test_target_plink Plink binary files for the target test population.
#' @param out_prefix Prefix for exported files. Recommended when plink files
#' are divided by chromosome or chunks, i.e., "chr1" or "chr2", and CTSLEB
#' script is dispatched in parallel.
#' @param r2_vec r square vector. Used by Plink1.9 "--clump r^2" parameter.
#' Defaults to c(0.01,0.05,0.1,0.2,0.5,0.8).
#' @param wc_base_vec Base clumping window size. Defaults to c(50,100)
#' @param mem Plink1.9 "--memory" parameter. Defaults to 8000 mb
#' @param threads Plink1.9 "--threads" parameter. Defaults to 2
#' @param params_farm List of plink parameters produced from SetParamsFarm()
#' @keywords plink1.9 clump
#' @export
#' @examples
#' data <- "data/"
#' results <- "test/"
#' plink19_exec <- "~/Apps/plink_v1.9/plink"
#' plink2_exec <- "~/Apps/plink2a/plink2"
#'
#' sum_EUR <- fread(paste0(data,"EUR_sumdata.txt"),header=T)
#' sum_AFR <- fread(paste0(data,"AFR_sumdata.txt"),header=T)
#'
#' Eur_ref_plinkfile <- paste0(data,"EUR_ref_chr22")
#' Afr_ref_plinkfile <- paste0(data,"AFR_ref_chr22")
#' Afr_test_plinkfile <- paste0(data,"AFR_test_chr22")
#' outprefix <- "chr22"
#'
#' PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
#'                           plink2_exec = plink2_exec)
#' prs_mat <- dimCT(results_dir = results,
#'                  sum_target = sum_AFR,
#'                  sum_ref = sum_EUR,
#'                  ref_plink = Eur_ref_plinkfile,
#'                  target_plink = Afr_ref_plinkfile,
#'                  test_target_plink = Afr_test_plinkfile,
#'                  out_prefix = outprefix,
#'                  params_farm = PRS_farm)

dimCT <- function(plink19_exec = 'plink',
                  plink2_exec = 'plink',
                  results_dir = "./",
                  sum_target,
                  sum_ref,
                  ref_plink,
                  target_plink,
                  test_target_plink,
                  out_prefix = as.null(),
                  r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
                  wc_base_vec = c(50,100),
                  memory= 8000,
                  threads = 2,
                  params_farm = as.null()) {

  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
    plink19_exec <- as.character(unlist(params_farm["plink19_exec"]))
    plink2_exec <-  as.character(unlist(params_farm["plink2_exec"]))
    r2_vec <- as.numeric(unlist(params_farm["r2_vec"]))
    wc_base_vec <- as.integer(unlist(params_farm["wc_base_vec"]))
    memory <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }

  sum_com <- AlignSum(sum_target = sum_AFR,
                      sum_ref = sum_EUR)
  assign("sum_com", sum_com, envir = .GlobalEnv)
  write_list <- SplitSum(x = sum_com,
                         results_dir = results_dir,
                         write_tables = TRUE)
  assign("write_list", write_list, envir = .GlobalEnv)

  ref_splitfile <- unlist(write_list["ref_split_file"])
  target_splitfile <- unlist(write_list["target_split_file"])
  snp_list <- RunClump(params_farm = params_farm,
                       plink19_exec = plink19_exec,
                       ref_plink = ref_plink,
                       target_plink = target_plink,
                       ref_splitfile = ref_splitfile,
                       target_splitfile = target_splitfile,
                       out_prefix = out_prefix,
                       results_dir = results_dir)
  assign("snp_list", snp_list, envir = .GlobalEnv)


  plink_list <- PreparePlinkFile(params_farm = params_farm,
                                  snp_list = snp_list,
                                  sum_com = sum_com,
                                  results_dir = results_dir)

  file_list <- helper_PreparePlinkFile(plink_list = plink_list,
                                       results_dir = results_dir)
  plink_list <- c(plink_list,file_list)
  assign("plink_list", plink_list, envir = .GlobalEnv)

  prs_mat <- PRSscore(params_farm = params_farm,
                      plink2_exec = plink2_exec,
                      bfile = test_target_plink,
                      plink_list = plink_list,
                      threads = threads,
                      memory = memory,
                      out_prefix = out_prefix,
                      results_dir = results_dir)
  assign("prs_mat", prs_mat, envir = .GlobalEnv)
  return(prs_mat)
}
