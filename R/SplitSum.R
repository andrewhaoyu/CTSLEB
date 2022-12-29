#' Inputs an AlignSum() object and splits SNPs into two groups
#' @description
#' Split the SNPs into two groups: 1. SNPs with p-values in the ref population
#' smaller than the target population. 2. SNPs with p-values in the target
#' population smaller than the ref population.
#'
#' @param x sum_com data.frame object produced by AlignSum().
#' @param results_dir
#' @param write_tables
#' @return  A list with the two data.frame objects of SNPs with p-values in the
#' ref population smaller than the target population in the [[1]] slot.
#' SNPs with p-values in the target population smaller than the ref population
#' in the [[2]] slot.
#' @export
#' @examples
#' sum_com <- AlignSum(sum_tar = sum_AFR, sum_other = sum_EUR)
#' split_list <- SplitSum(sum_com)

SplitSum <- function(x,
                     results_dir,
                     write_tables = TRUE){
  sum.com <- x
  sum.com_select <- sum.com %>%
    mutate(split_ind =
             ifelse(
               (P<P_ref)|is.na(P_ref),T,F)
    )%>%
    select(SNP,P,P_ref,split_ind)

  sum.com_select_other_ref <- sum.com_select %>%
    filter(split_ind==F) %>%
    select(SNP,P_ref) %>%
    rename(P = P_ref)

  sum.com_select_tar_ref <- sum.com_select %>%
    filter(split_ind==T) %>%
    select(SNP,P)

  sum_list <- list(sum.com_select_other_ref,
                   sum.com_select_tar_ref)
  if (write_tables) {
    write_list <- WriteSplitTables(x = sum_list,
                             results_dir = results_dir,
                             ref_split_file = "sum_ref.txt",
                             target_split_file = "sum_target.txt")

    return(write_list)
  } else {
    print(paste0("WriteSplitTables() was not performed"))
    return(sum_list)
  }

  return(sum_list)

}
