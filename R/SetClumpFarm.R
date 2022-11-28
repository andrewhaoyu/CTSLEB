#' set Plink clumping parameters
#' @description
#' This function sets the plink1.9 clumping parameters
#' flags values, r square and base window size vectors
#' @param plink19exec path and name to plink1.9 executable. Default "plink"
#' @param wc_base_vec Base window size. Defaults to c(50,100)
#' @param r2_vec r square vector. Defaults to c(0.01,0.05,0.1,0.2,0.5,0.8)
#' @param threads --threads for plink1.9 Defaults to 2
#' @param mem --memory for plink1.9 Defaults to 8000 mb
#' @param ref_out folder and ref file prefix for plink output
#' @param target_out folder and target file prefix for plink output
#' @return A named list
#' @keywords plink1.9 clump
#' @export
#' @examples

#' clump_farm <- SetClumpFarm(plink19exec = /Apps/plink1.9/plink,
#'                            wc_base_vec = c(50,100),
#'                            r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
#'                            threads = 2,
#'                            mem = 8000
#' )


SetClumpFarm <- function(plink19exec = 'plink',
                         wc_base_vec=c(50,100),
                         r2_vec=c(0.01,0.05,0.1,0.2,0.5,0.8),
                         threads = 2,
                         mem = 8000)
  {
  # clump_farm_ld_names <- c("plink19exec",
  #                          "ref_plink",
  #                          "target_plink",
  #                          "ref_split_file",
  #                          "target_split_file",
  #                          "r2_vec",
  #                          "wc_base_vec",
  #                          "threads",
  #                          "mem",
  #                          "ref_out",
  #                          "target_out")
  # clump_farm_ld_values <- list(plink19exec,
  #                              ref_plink,
  #                              target_plink,
  #                              ref_split_file,
  #                              target_split_file,
  #                              r2_vec,
  #                              wc_base_vec,
  #                              threads,
  #                              mem,
  #                              ref_out,
  #                              target_out)
  clump_farm_ld_names <- c("plink19exec",
                           "r2_vec",
                           "wc_base_vec",
                           "threads",
                           "mem")
  clump_farm_ld_values <- list(plink19exec,
                               r2_vec,
                               wc_base_vec,
                               threads,
                               mem)
  clump_farm <- setNames(clump_farm_ld_values, clump_farm_ld_names)
  print(clump_farm)
  return(clump_farm)
}
