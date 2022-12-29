#' Creates the two data tables created by SplitSum()
#' @description
#' Creates two data.frame objects from SplitSum() result and writes them to
#' files. These files are inputs to the Plink 1.9 --clump flag for LD clumping.
#' @param x List containing data.frame objects created by SplitSum()
#' @param refFile export path and prefix file name for ref data.frame
#' Output file will be named "<refFile>_ref"
#' @param targetFile export path and prefix file name for target data.frame
#' Output file will be named "<targetFile>_target"
#' @return The ref and target aligned data.frames produced by SplitSum() are
#' written to default files named "sum_ref.txt" and "sum_target.txt" respectively
#' @examples
#' sum_com_split <- SplitSum(sum_com)
#' WriteSplitTables(x=sum_com_split, refFile = "sum", targetFile = "sum")

WriteSplitTables <- function(x,
                             results_dir,
                             ref_split_file = "sum_ref.txt",
                             target_split_file = "sum_target.txt")
  {
  print("executing WriteSplitTables()... ")
  ref_split_file <- paste0(results_dir, ref_split_file)
  target_split_file <- paste0(results_dir,target_split_file)
  print("creating sum_ref object")
  sum_ref <- x[[1]]
  #assign("sum_ref", x[[1]], envir = .GlobalEnv)
  print("creating sum_target object")
  #assign("sum_target", x[[2]], envir = .GlobalEnv)
  sum_target <- x[[2]]
  assign("ref_split_file", ref_split_file, envir = .GlobalEnv)
  assign("target_split_file", target_split_file, envir = .GlobalEnv)
  write.table(sum_ref,ref_split_file,col.names = T,row.names = F,quote=F)
  write.table(sum_target,target_split_file,col.names = T,row.names = F,quote=F)
  names <- c("ref_split_file",
             "target_split_file")
  values <- list(ref_split_file,
                 target_split_file)
  write_list <- setNames(values, names)
  print("WriteSplitTables() complete ... ")
  return(write_list)
}
