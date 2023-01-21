#' Inputs an AlignSum() object and splits SNPs into two groups
#' @description
#' Split the SNPs into two groups: 1. SNPs with p-values in the ref population
#' smaller than the target population. 2. SNPs with p-values in the target
#' population smaller than the ref population.
#' @param x sum_com Data.frame object produced by AlignSum().
#' @param results_dir Folder or directory to export files. Mandatory when
#' 'write_tables = TRUE'. Default file names are sum_ref.txt and sum_target.txt.
#' @param write_tables Writes two files to the 'results_dir' used by the plink1.9
#' '--clump' parameter. Default is TRUE.
#' @return  A list with the two data.frame objects.SNPs with p-values in the
#' ref population smaller than the target population in the [[1]] slot.
#' SNPs with p-values in the target population smaller than the ref population
#' in the [[2]] slot. If "write_tables = TRUE" the tables will be written as files
#' in the results_dir folder AND the global variables ref_split_file and
#' target_split_file will be created for the file locations.
#' @export
#' @examples
#' data <- "data/"
#' results <- "test/"
#' sum_EUR <- fread(paste0(data, "EUR_sumdata.txt"), header=T)
#' sum_AFR <- fread(paste0(data, "AFR_sumdata.txt"), header=T)
#' sum_com <- AlignSum(sum_tar = sum_AFR, sum_other = sum_EUR)
#' split_list <- SplitSum(sum_com = sum_com,
#'                        results_dir = results,
#'                        write_tables = TRUE)

SplitSum <- function(x,
                     results_dir,
                     write_tables = TRUE){
  print("executing SplitSum() ...")
  sum_com <- x
  sum_com_select <- sum_com %>%
    mutate(split_ind =
             ifelse(
               (P<P_ref)|is.na(P_ref),T,F)
    )%>%
    select(SNP,P,P_ref,split_ind)

  sum_com_select_other_ref <- sum_com_select %>%
    filter(split_ind==F) %>%
    select(SNP,P_ref) %>%
    rename(P = P_ref)

  sum_com_select_tar_ref <- sum_com_select %>%
    filter(split_ind==T) %>%
    select(SNP,P)

  sum_list <- list(sum_com_select_other_ref,
                   sum_com_select_tar_ref)

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
  print("SplitSum() complete ... ")
  return(sum_list)

}
