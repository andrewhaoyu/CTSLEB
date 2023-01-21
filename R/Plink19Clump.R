#' Plink19Clump
#'
#' This function sets the standard options LD clumping with plink1.9.
#' If additional options are required, please use the source() command
#' as demonstrated in the vignette.
#' @param plink19_exec Plink binary execute file. Default "plink".
#' @param clump PLINK-format association report(s) (text files with a header
#' line, a column containing variant IDs, and another column
#' containing p-values).
#' @param bfile Prefix for plink binary set (prefix.bed, prefix.bim, prefix.fam).
#' @param clump_p1 Significance threshold for a SNP to be included as an index SNP.
#' Default = 1 is used to include all SNPs for clumping.
#' @param clump_r2 LD r^2 threshold for clumping. SNPs with a higher r^2
#' with the index SNPs will be removed.
#' @param clump_kb Physical distance threshold for clumping. SNPs
#' within the radius of the index SNP are considered for clumping.
#' For PRS, clump_kb = base_window_size/clumping_r_square so that lower
#' clump_r2 thresholds have a larger clumping window size.
#' @param threads Maximum number of concurrent threads.
#' @param memory Primary workspace memory.
#' @param out Output folder and file prefix for Plink1.9. Output files will end
#' with .clumped
#' @keywords plink1.9 clump
#' @usage Plink19Clump(plink19_exec, clump, bfile, out, params_farm)
#' @export
#' @usage Plink19Clump(plink19_exec,bfile,clump,clump_p1,clump_r2,clump_kb,
#' threads,memory,out,params_farm)


Plink19Clump <- function(plink19_exec = "plink",
                         bfile,
                         clump,
                         clump_p1 = 1,
                         clump_r2,
                         clump_kb,
                         threads = 4,
                         memory = 8000,
                         out = "./",
                         params_farm=as.null()){

  if (is.null(params_farm)) {
    print("Plink19Clump() no params_farm")
  } else {
    print("Plink19Clump() params_farm list will be used")
    mem <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
  }
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
