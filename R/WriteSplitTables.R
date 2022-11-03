#' Creates the two data tables created by SplitSum() and writes them to
#' files used by Plink 1.9 for LD clumping.
#'
#' @param x List containing data.tables created by SplitSum()
#' @param refFile export path and prefix file name for ref data.table.
#' Output file will be named "<refFile>_ref"
#' @param targetFile export path and prefix file name for target data.table.
#' Output file will be named "<targetFile>_target"


WriteSplitTables <- function(x,
                             refFile = "sum",
                             targetFile = "sum") {
  print(refFile)
  print(targetFile)
}
