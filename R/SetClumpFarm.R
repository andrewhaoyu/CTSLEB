#' set Plink clumping parameters
#' @description
#' This function sets the plink1.9 clumping parameters
#' flags values, r square and base window size vectors
#' @param plink19exec path and name to plink1.9 executable. Default "plink"
#' @param ref_plink Location and prefix for reference plink binary data set
#' @param target_plink Location and prefix for target plink binary data set
#' @param ref_split Reference snps data table generated from SplitSum().
#' Defaults to sum_other_ref.
#' @param target_split  Target snps data table generated from SplitSum().
#' Defaults to sum_tar_ref
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
#' temp.dir <- test/temp.dir
#' system("mkdir -p test/temp.dir")
#' clump_farm <- SetClumpFarm(plink19exec = /Apps/plink1.9/plink,
#'                            ref_plink = "temp/plink/Eur",
#'                            target_plink = "temp/plink/Afr",
#'                            ref_split = "temp/Eur_split",
#'                            target_split = "temp/Afr_split",
#'                            wc_base_vec = c(50,100),
#'                            r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
#'                            threads = 2,
#'                            mem = 8000,
#'                            ref_out = 'temp/ref'
#'                            target_out = 'temp/target'
#' )


SetClumpFarm <- function(plink19exec = 'plink',
                         ref_plink,
                         target_plink,
                         ref_split,
                         target_split,
                         wc_base_vec=c(50,100),
                         r2_vec=c(0.01,0.05,0.1,0.2,0.5,0.8),
                         threads = 2,
                         mem = 8000){
  clump_farm_ld_names <- c("plink19exec",
                           "ref_plink",
                           "target_plink",
                           "ref_split",
                           "target_split",
                           "r2_vec",
                           "wc_base_vec",
                           "threads",
                           "mem",
                           "ref_out",
                           "target_out")
  clump_farm_ld_values <- list(plink19exec,
                               ref_plink,
                               target_plink,
                               ref_split,
                               target_split,
                               r2_vec,
                               wc_base_vec,
                               threads,
                               mem)
  clump_farm <- setNames(clump_farm_ld_values, clump_farm_ld_names)
  print(clump_farm)
  return(clump_farm)
}
