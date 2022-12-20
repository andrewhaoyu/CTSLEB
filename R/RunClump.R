#'
#' Perform clumping with plink1.9
#' @param plink19_exec Plink binary execute file. Default "plink".
#' @param sum_ref
#' @param sum_target
#' @param results_dir
#' @param ref_plink
#' @param target_plink
#' @param ref_outprefix
#' @param target_outprefix
#' @param r2_vec
#' @param wc_base_vec
#' @param mem
#' @param threads
#' @param params_farm

#' @keywords plink1.9 clump
#' @usage dimCT(plink19_exec, params_farm)

RunClump <- function(plink19_exec=plink19_exec,
                     ref_plink = ref_plink,
                     target_plink = target_plink,
                     ref_outprefix = "ref",
                     ref_outprefix = "target",
                     results_dir = results_dir,
                     ref_split_file = ref_split_file,
                     target_split_file = target_split_file,
                     r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
                     wc_base_vec = c(50,100),
                     mem = 8000,
                     threads = 2,
                     params_farm = as.null()) {
  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
    plink19_exec <- as.character(unlist(params_farm["plink19_exec"]))
    r2_vec <- as.vector(unlist(params_farm["r2_vec"]))
    wc_base_vec <- as.vector(unlist(params_farm["wc_base_vec"]))
    mem <- as.character(unlist(params_farm["mem"]))
    threads <- as.character(unlist(params_farm["threads"]))
  }

  snp_list <-list()
  temp <- 1
  temp.dir <- paste0(results_dir,"temp/")
  for(r_ind in 1:length(r2_vec)){
    wc_vec <- wc_base_vec/r2_vec[r_ind]
    for(w_ind in 1:length(wc_vec)){
      pthr <-1
      r2thr <- r2_vec[r_ind]
      kbpthr <- wc_vec[w_ind]
      ref_outfile <- paste0(temp.dir,ref_outprefix, "_CT_rind_",r_ind,"_wcind_",w_ind)
      target_outfile <- paste0(temp.dir,target_outprefix, "_CT_rind_",r_ind,"_wcind_",w_ind)
      Plink19Clump(plink19_exec = plink19_exec,
                   bfile = ref_plink,
                   clump = ref_split_file,
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
                         clump = target_split_file,
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
      print(paste0("binding rows for ",
                   ref_outfile,
                   ".clumped and ",
                   target_outfile,
                   ".clumped"))
      LD  <- rbind(LD_ref,LD_target)
      snp_list[[temp]] <- LD
      print(paste0("creating snp list for clump_r2_",
                   r2thr,
                   "_ws_",
                   kbpthr))
      names(snp_list[[temp]]) <- paste0("clump_r2_",
                                        r2thr,
                                        "_ws_",
                                        kbpthr)
      temp <- temp + 1
    }
  }
  assign("snp_list", snp_list, envir = .GlobalEnv)
}
