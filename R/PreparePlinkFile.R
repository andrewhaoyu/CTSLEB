#' Prepare the data frames for PLINK2 --score routine
#' @description
#' Prepares the scores, p_values, unique_infor and q_range data frames required
#' for PLINK2 --score routine
#' @param snp_list The snp_list object from RunClump().
#' @param sum_com the sum_com result from AlignSum().
#' @param pthres vector of p-value thresholds. Default
#' c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
#' @param results_dir Folder or directory to export files.
#' @return List with data frames scores [[1]], p_values [[2]],
#' unique_infor [[3]] and q_range [[4]]. The scores data.frame first column
#' contains the unique SNPs for the LD clumping combinations. The second column
#' is the effect allele. The third through last columns contain the regression
#' coefficients of the target population for SNPs after the LD-clumping combination.
#' The coefficients is set to 0 if a SNP doesn't exist in that clumping combination.
#' The p_value data.frame first column is the same as scores. The second column
#' contains the p-values for the SNPs from the target population GWAS. The
#' unique_infor data.frame contains the information for SNPs after the clumping
#' step; SNP, CHR, BP, A1 (effect_allele), (BETA, SE, P) for the target population,
#' (BETA_other, SE_other,P_other) for the reference  population
#' @export
#' @examples
#' data <- "data/"
#' results <- "test/"
#' temp <- paste0(results,"temp/")
#' plink19_exec <- "~/Apps/plink_v1.9/plink"
#' plink2_exec <- "~/Apps/plink2a/plink2"
#'
#' sum_EUR <- fread(paste0(data,"EUR_sumdata.txt"),header=T)
#' sum_AFR <- fread(paste0(data,"AFR_sumdata.txt"),header=T)
#'
#' Eur_ref_plinkfile <- paste0(data,"EUR_ref_chr22")
#' Afr_ref_plinkfile <- paste0(data,"AFR_ref_chr22")
#' Afr_test_plinkfile <- paste0(data,"AFR_test_chr22")
#' outprefix <- "chr22"
#'
#' PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
#'                           plink2_exec = plink2_exec)
#' sum_com <- AlignSum(sum_target = sum_AFR,
#'                     sum_ref = sum_EUR)
#' write_list <- SplitSum(x = sum_com,
#'                        results_dir = results)
#' snp_list <- RunClump(params_farm = PRS_farm,
#'                      plink19_exec = plink19_exec,
#'                      ref_plink = Eur_ref_plinkfile,
#'                      target_plink = Afr_ref_plinkfile,
#'                      ref_splitfile = ref_split_file,
#'                      target_splitfile = target_split_file,
#'                      out_prefix = outprefix,
#'                      results_dir = results)
#' plink_list <- PreparePlinkFile(params_farm = PRS_farm,
#'                                snp_list = snp_list,
#'                                sum_com = sum_com,
#'                                results_dir = results_dir)

PreparePlinkFile <- function(snp_list = snp_list,
                             sum_com,
                             pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                             results_dir,
                             params_farm=as.null())
  {
  print("executing PreparePlinkFile()... ")
  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }

  #create unique SNP list by combind LD clumping results under different parameters

  unique_id <- unique(rbindlist(snp_list,use.name =FALSE))
  names(unique_id) <- "SNP"

  #align the regression coefficients for these SNPs from the sum stat

  unique_infor <- left_join(unique_id,sum_com,by="SNP")
  print("unique_infor dataframe complete")

  #create a coefficient matrix

  n_col <- length(snp_list)
  n_row <- nrow(unique_infor)
  beta_mat <- matrix(unique_infor$BETA,nrow =n_row,ncol =n_col)
  names <- rep("c",n_col)
  temp <- 1
  for(ldx in seq(n_col)){
    LD <- snp_list[[ldx]]
    names(LD) <- "SNP"
    idx <- which(unique_infor$SNP%in%LD$SNP==F)
    beta_mat[idx,ldx] <- 0
    names[ldx] <- names(snp_list[[ldx]])
  }
  colnames(beta_mat) <- names
  scores <- data.frame(SNP = unique_id,A1 = unique_infor$A1,beta_mat)
  print("scores dataframe complete")
  p_values <- data.frame(SNP = unique_id,P = unique_infor$P)
  print("p_values dataframe complete")
  q_range <- CreateQRange(pthres)

  names <- c("scores",
             "p_values",
             "unique_infor",
             "q_range")
  values <- list(scores,
                 p_values,
                 unique_infor,
                 q_range)
  list <- setNames(values, names)
  print("PreparePlinkFile() complete ...")

  return(list)

}
