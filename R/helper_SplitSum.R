#' Helper function for AlignSum()
#' @description
#' SplitSum = TRUE, call the SplitSum() and WriteSplitTables() functions and
#' create the sum_com global variable object.
#' SplitSum = FALSE, only create the sum_com global variable object.
#' Not exported
#' @para x sum_com Object from AlignSum().
#' @param results_dir
#' @return SplitSum() and WriteSplitTables() objects and files

helper_SplitSum <- function(sum_com, results_dir = results_dir)
  # ref_split_file = ref_split_file
  # target_split_file = target_split_file
  {
  split_list <- SplitSum(sum_com)
  WriteSplitTables(x = split_list, results_dir = results_dir)
  # ref_split_file = ref_split_file,
  # target_split_file = target_split_file)
}
#@param ref_split_file If SplitSum = TRUE, this is the file path and name for
#the reference snps produced from the SplitSum() function.  The default name
#"sum_ref.txt" is passed to WriteSplitTables().
#target_split_file If SplitSum = TRUE, this is the file path and name
#for the target snps produced from the SplitSum() function. The default
#name "sum_target.txt" is passed to WriteSplitTables().
