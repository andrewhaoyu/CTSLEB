#' set Plink clumping parameters
#' @description
#' This function sets the plink1.9 clumping and plink2a thresholding parameters
#' and flags
#' @param plink19_exec path and name to plink1.9 executable. Default "plink"
#' @param plink2_exec path and name to plink2 executable. Default "plink2"
#' @param wc_base_vec Base window size. SetClumpFarm will automaticaly create the
#' global 'wc_base_vec' since it is required in other workflow sections. Defaults
#' to c(50,100)
#' @param r2_vec r square vector. SetClumpFarm will automatically create the
#' global 'r2_vec' since it is required in other workflow sections. Defaults to
#' c(0.01,0.05,0.1,0.2,0.5,0.8).
#' @param pthres  description
#' @param threads --threads for plink1.9 Defaults to 2
#' @param mem --memory for plink1.9 Defaults to 8000 mb
#' @param ref_out folder and ref file prefix for plink output
#' @param target_out folder and target file prefix for plink output
#' @return A named list
#' @keywords plink1.9 clump
#' @export
#' @usage SetParamsFarm(plink19exec, wc_base_vec, r2_vec, pthres, threads,
#' mem)

SetParamsFarm <- function(plink19_exec = 'plink',
                       plink2_exec = 'plink2',
                       wc_base_vec = c(50,100),
                       r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
                       pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                       threads = 2,
                       mem = 8000)
  {

  names <- c("plink19_exec",
                      "plink2_exec",
                         "r2_vec",
                         "wc_base_vec",
                         "pthres",
                         "threads",
                         "mem")
  values <- list(plink19_exec,
                          plink2_exec,
                          r2_vec,
                          wc_base_vec,
                          pthres,
                          threads,
                          mem)
  params_farm <- setNames(values, names)
  assign("r2_vec", r2_vec, envir = .GlobalEnv)
  assign("wc_base_vec", wc_base_vec, envir = .GlobalEnv)
  assign("pthres", pthres, envir = .GlobalEnv)
  print(params_farm)
  return(params_farm)
}
