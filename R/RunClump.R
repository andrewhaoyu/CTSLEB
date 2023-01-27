#'
#' LD clumping with plink1.9
#' @param plink19_exec Plink1.9  binary execute file. Default "plink".
#' @param results_dir Folder or directory to export files.
#' @param ref_plink Plink binary files for the reference population.
#' @param target_plink Binary plink files for the target population.
#' @param ref_splitfile File name of reference file from SplitSum().
#' @param target_splitfile File name of target file from SplitSum().
#' @param out_prefix Prefix for exported files. Recommended when plink files
#' are divided by chromosome or chunks, i.e., "chr1" or "chr2", and CTSLEB
#' script is dispatched in parallel.
#' @param r2_vec r square vector. Used by Plink1.9 "--clump r^2" parameter.
#' Defaults to c(0.01,0.05,0.1,0.2,0.5,0.8).
#' @param wc_base_vec Base clumping window size. Defaults to c(50,100)
#' @param mem Plink1.9 "--memory" parameter. Defaults to 8000 mb
#' @param threads Plink1.9 "--threads" parameter. Defaults to 2
#' @param params_farm List of plink parameters produced from SetParamsFarm()
#' @export
#' @keywords plink1.9 clump
#' @examples
#' data <- "data/"
#' results <- "test/"
#' temp <- paste0(results,"temp/")
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
#' sum_com <- AlignSum(sum_target = sum_AFR,
#'                     sum_ref = sum_EUR)
#' write_list <- SplitSum(x = sum_com,
#'                        results_dir = results)
#' snp_list <- RunClump(params_farm = PRS_farm,
#'                      plink19_exec = plink19_exec,
#'                      ref_plink = Eur_ref_plinkfile,
#'                      target_plink = Afr_ref_plinkfile,
#'                      ref_splitfile = ref_split_file,
#'                      target_splitfile = target_split_file,
#'                      out_prefix = outprefix,
#'                      results_dir = results)

RunClump <- function(plink19_exec,
                     ref_plink,
                     target_plink,
                     out_prefix = as.null(),
                     results_dir,
                     ref_splitfile,
                     target_splitfile,
                     r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
                     wc_base_vec = c(50,100),
                     mem = 8000,
                     threads = 2,
                     params_farm = as.null()) {
  print("executing RunClump()... ")
  if (is.null(params_farm)) {
    print("RunClump() no params_farm")
  } else {
    print("RunClump() params_farm list will be used")
    plink19_exec <- as.character(unlist(params_farm["plink19_exec"]))
    r2_vec <- as.numeric(unlist(params_farm["r2_vec"]))
    wc_base_vec <- as.integer(unlist(params_farm["wc_base_vec"]))
    mem <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
  }

  temp.dir <- paste0(results_dir,"temp/")
  if (is.null(out_prefix)) {
    print("not out_prefix")
    out.prefix <- ""
  } else {
    print(paste0("out_prefix: ", out_prefix))
    out.prefix <- paste0(out_prefix, "_")
  }

  snp_list <-list()
  temp <- 1

  for(r_ind in 1:length(r2_vec)){
    wc_vec <- wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      pthr <-1
      r2thr <- r2_vec[r_ind]
      kbpthr <- wc_vec[w_ind]
      ref_outfile <- paste0(temp.dir, out.prefix, "ref_CT_rind_",r_ind,"_wcind_",w_ind)
      target_outfile <- paste0(temp.dir, out.prefix,"target_CT_rind_",r_ind,"_wcind_",w_ind)
      Plink19Clump(plink19_exec = plink19_exec,
                   bfile = ref_plink,
                   clump = ref_splitfile,
                   clump_p1 = pthr,
                   clump_r2 = r2thr,
                   clump_kb = kbpthr,
                   threads = threads,
                   memory = mem,
                   out = ref_outfile,
                   params_farm = params_farm
      )
      Plink19Clump(plink19_exec = plink19_exec,
                   bfile = target_plink,
                   clump = target_splitfile,
                   clump_p1 = pthr,
                   clump_r2 = r2thr,
                   clump_kb = kbpthr,
                   threads = threads,
                   memory = mem,
                   out = target_outfile,
                   params_farm = params_farm
      )
      LD_ref <-fread(paste0(ref_outfile,
                            ".clumped"))[,3,drop=F]
      LD_target <-fread(paste0(target_outfile,
                               ".clumped"))[,3,drop=F]
      # print(paste0("binding rows for ",
      #              ref_outfile,
      #              ".clumped and ",
      #              target_outfile,
      #              ".clumped"))
      LD  <- rbind(LD_ref,LD_target)
      snp_list[[temp]] <- LD
      # print(paste0("creating snp list for clump_r2_",
      #              r2thr,
      #              "_ws_",
      #              kbpthr))
      names(snp_list[[temp]]) <- paste0("clump_r2_",
                                        r2thr,
                                        "_ws_",
                                        kbpthr)
      temp <- temp + 1
    }
  }
  print("RunClump() complete ...")
  return(snp_list)
}
