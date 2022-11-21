#' plink19Clump
#'
#' This function sets the standard options LD clumping with plink1.9.
#' If additional options are required, please use the source() command
#' as demonstrated in the vignette.
#' @param plink19_exec Plink binary execute file. Default "plink".
#' @param clump PLINK-format association report(s) (text files with a header
#' line, a column containing variant IDs, and another column
#' containing p-values)
#' @param bfile prefix for plink binary set (prefix.bed, prefix.bim, prefix.fam)
#' @param clump_p1 Significance threshold for a SNP to be included as an index SNP.
#' Default = 1 is used to include all SNPs for clumping
#' @param clump_r2 LD r^2 threshold for clumping. SNPs with a higher r^2
#' with the index SNPs will be removed.
#' @param clump_kb physical distance threshold for clumping. SNPs
#' within the radius of the index SNP are considered for clumping.
#' For PRS, clump_kb = base_window_size/clumping_r_square so that lower
#' clump_r2 thresholds have a larger clumping window size
#' @param threads maximum number of concurrent threads
#' @param memory primary workspace memory
#' @param out out file prefix. Output files end with .clumped
#' @keywords plink1.9 clump
#' @usage Plink19Clump(plink19_exec, clump_farm)
#' @export
#' @examples
#' r2_vec <- c(0.01,0.05,0.1,0.2,0.5,0.8)
#' wc_base_vec <- c(50,100)
#' bfile <- "data/EUR_ref_chr22"
#' clump <- "temp/sum_other_ref"
#' clump_r2 <- r2_vec[1]
#' clump_kb <- wc_base_vec[1]/clump_r2
#' plink19_clump(plink19_exec,
#' bfile = bfile,
#' clump = clump,
#' clump_r2 = clump_r2,
#' clump_kb = clump_kb,
#' out = temp/ref)

Plink19Clump <- function(plink19_exec = "plink",
                         bfile,
                         clump,
                         clump_p1 = 1,
                         clump_r2,
                         clump_kb,
                         threads = 4,
                         memory = 8000,
                         out){
  system(paste0(plink19_exec, " ",
                "--bfile ", bfile, " ",
                "--clump ", clump, " ",
                "--clump-p1 ", clump_p1, " ",
                "--clump-r2 ", clump_r2, " ",
                "--clump-kb ", clump_kb, " ",
                "--threads ", threads, " ",
                "--memory ", memory, " ",
                "--out ", out)
  )
}
