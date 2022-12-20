#' Plink19Clump
#'
#' This function sets the standard options LD clumping with plink1.9.
#' If additional options are required, please use the source() command
#' as demonstrated in the vignette.
#' @param plink19_exec Plink binary execute file. Default "plink".
#' @param clump PLINK-format association report(s) (text files with a header
#' line, a column containing variant IDs, and another column
#' containing p-values)
#' @param bfile prefix for plink binary set (prefix.bed, prefix.bim, prefix.fam)
#' @param clump_r2 LD r^2 threshold for clumping. SNPs with a higher r^2
#' with the index SNPs will be removed.
#' @param clump_kb physical distance threshold for clumping. SNPs
#' within the radius of the index SNP are considered for clumping.
#' For PRS, clump_kb = base_window_size/clumping_r_square so that lower
#' clump_r2 thresholds have a larger clumping window size
#' @param threads maximum number of concurrent threads
#' @param memory primary workspace memory
#' @param out output folder and file prefix for Plink1.9. Output files will end
#' with .clumped
#' @keywords plink1.9 clump
#' @usage Plink19Clump(plink19_exec, clump, bfile, out, params_farm)
#' @export
#' @examples
#' r2_vec <- c(0.01,0.05,0.1,0.2,0.5,0.8)
#' wc_base_vec <- c(50,100)
#' bfile <- "data/EUR_ref_chr22"
#' clump <- "temp/sum_ref.txt"
#' clump_r2 <- r2_vec[1]
#' clump_kb <- wc_base_vec[1]/clump_r2
#' clump_out <- paste0(temp.dir, "ref")
#' plink19_clump(plink19_exec,
#'               bfile = refFile,
#'               clump = clump,
#'               clump_r2 = clump_r2,
#'               clump_kb = clump_kb,
#'               out = clump_out)

PRSclump <- function(plink19_exec = "plink",
                     ref_plink,
                     target_plink,
                     results_dir,
                     ref_split_file,
                     target_split_file,
                     ref_clump_out,
                     target_clump_out,
                     r2_vec = c(0.01,0.05,0.1,0.2,0.5,0.8),
                     wc_base_vec = c(50,100),
                     threads = 4,
                     memory = 8000,
                     params_farm=as.null()){

  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
  }

  helper_clump(ref_plink = ref_plink,
               target_plink = target_plink,
               results_dir = results_dir,
               ref_split_file = ref_split_file,
               target_split_file = target_split_file,
               ref_clump_out = ref_clump_out,
               target_clump_out = target_clump_out,
               r2_vec = r2_vec,
               wc_base_vec = c(50,100),
               threads = 4,
               memory = 8000,
               params_farm = params_farm)
}
