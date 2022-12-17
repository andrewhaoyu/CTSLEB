#' Helper function for AlignSum()
#' @description
#' SplitSum = TRUE, call the SplitSum() and WriteSplitTables() functions and
#' create the sum_com global variable object.
#' SplitSum = FALSE, only create the sum_com global variable object.
#' Not exported
#' @param x TRUE default
#' @para sum_com Object from AlignSum().
#' @param ref_split_file If SplitSum = TRUE, this is the file path and name for
#' the reference snps produced from the SplitSum() function.  The default name
#' "sum_ref.txt" is passed to WriteSplitTables().
#' @param target_split_file If SplitSum = TRUE, this is the file path and name
#' for the target snps produced from the SplitSum() function. The default
#' name "sum_target.txt" is passed to WriteSplitTables().
#' @return SplitSum() and WriteSplitTables() objects and files
#'
helper_SplitSum <- function(sum_com,
                            results_dir,
                            ref_split_file,
                            target_split_file
                            )
  {
    split_list <- SplitSum(sum_com)

    WriteSplitTables(x = split_list,
                     results_dir = results_dir,
                     ref_split_file = ref_split_file,
                     target_split_file = target_split_file)
}
