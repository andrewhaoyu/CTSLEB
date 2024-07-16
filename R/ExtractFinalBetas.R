#' Clean the PRS's for the Tune and Validation set for the Super Learner step.
#'
#' @param Tune_PRS Cleaned Tune PRS's that were used in the Super Learner step. Remove IID and FID, i.e. the first two columns.
#' @param Predicted_Tune_Y Raw predicted tune response from the Super Learner.
#' @param prs_mat_eb prs_mat_eb from the EB step.
#' @param unique_infor_post Global environment object from prior steps.
#' @param pthres pthres vector set in earlier steps as a parameter.
#'
#' @return Named list of Cleaned Tune PRS and Cleaned Validation PRS.
#' @export

ExtractFinalBetas <- function(Tune_PRS,Predicted_Tune_Y,prs_mat_eb,unique_infor_post,pthres){
  Tune_PRS <- as.matrix(Tune_PRS)
  Y <- matrix(Predicted_Tune_Y,ncol = 1)
  Beta_Star <- matrix(unname(coef(lm(Predicted_Tune_Y~Tune_PRS))),ncol = 1)
  Beta_Star <- Beta_Star[-1,,drop = FALSE]

  score_full <- matrix(NA,ncol = ncol(prs_mat_eb),nrow = nrow(unique_infor_post))
  names_score_full <- NULL
  count <- 1
  for(i in 1:length(pthres)){
    for(j in 1:length(pthres)){
      for(k in 1:24){
        tmp <- matrix(0,nrow = nrow(unique_infor_post),ncol = 1)
        tmp[which(unique_infor_post$P < pthres[j] | unique_infor_post$P_ref < pthres[i]),1] <- plink_list_eb$scores_eb[which(unique_infor_post$P < pthres[j] | unique_infor_post$P_ref < pthres[i]),k + 2]
        names_score_full <- c(names_score_full,paste0(colnames(plink_list_eb$scores_eb)[k + 2],"_p_other_",pthres[i],"_p_tar_",pthres[j]))
        score_full[,count] <- tmp[,1]
        count <- count + 1
      }
    }
  }

  score_full <- score_full[,colnames(prs_mat_eb)[-c(1,2)] %in% colnames(Tune_PRS)]
  names_eb <- colnames(prs_mat_eb)[colnames(prs_mat_eb) %in% colnames(Tune_PRS)]
  score_full <- score_full[,match(colnames(Tune_PRS),names_eb)]
  final <- data.frame(SNP = unique_infor_post$SNP,A1 = unique_infor_post$A1,BETA = score_full%*%Beta_Star)
  return(final)
}
