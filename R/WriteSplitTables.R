#' Creates the two data tables created by SplitSum()
#' @description
#' Creates the two data tables created by SplitSum() and writes them to
#' files used by Plink 1.9 for LD clumping.
#' @param x List containing data.tables created by SplitSum()
#' @param refFile export path and prefix file name for ref data.table.
#' Output file will be named "<refFile>_ref"
#' @param targetFile export path and prefix file name for target data.table.
#' Output file will be named "<targetFile>_target"
#' @return Two data.table objects named "sum_ref" for the reference snps and
#' "sum_target" for the target snps created by SplitSum().
#' @export
#' @examples
#' sum_com_split <- SplitSum(sum_com)
#' WriteSplitTables(x=split_com_split, refFile = "sum", targetFile = "sum")

WriteSplitTables <- function(x,
                             refFile = "sum",
                             targetFile = "sum") {
  ref_out <- paste0(tools::file_path_sans_ext(baseFile),"_ref")
  target_out <- paste0(tools::file_path_sans_ext(targetFile),"_target")
  print("creating sum_ref")
  assign("sum_ref", x[[1]], envir = .GlobalEnv)
  print("creating sum_target")
  assign("sum_target", x[[2]], envir = .GlobalEnv)

  write.table(sum_ref,paste0(ref_out),col.names = T,row.names = F,quote=F)
  write.table(sum_target,paste0(target_out),col.names = T,row.names = F,quote=F)
}
