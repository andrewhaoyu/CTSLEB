#' Clean the PRS's for the Tune and Validation set for the Super Learner step.
#'
#' @param Tune_PRS PRS's for the Tune Data built using the prs_mat_eb object. First column is IID and second column is FID.
#' @param Tune_Y The observed phenotypes for the Tune Data.
#' @param Validation_PRS PRS's for the Validation Data built using the prs_mat_eb object. First column is IID and second column is FID.
#' @param family Distribution of the response, either "Gaussian" or "Binary".
#'
#' @return Named list of Cleaned Tune PRS and Cleaned Validation PRS. First two columns of each data.frame are IID and FID.
#' @export

PRS_Clean <- function(Tune_PRS,Tune_Y,Validation_PRS,family = "Gaussian"){

  if(sum(family %in% c("Gaussian","Binary") != 1)){
    stop("family has to be one of Gaussian or Binary")
  }

  corcutoff <- 0.98

  Tune_IID <- Tune_PRS[,1]
  Tune_FID <- Tune_PRS[,2]
  Tune_PRS <- Tune_PRS[,-c(1,2)]

  Validation_IID <- Validation_PRS[,1]
  Validation_FID <- Validation_PRS[,2]
  Validation_PRS <- Validation_PRS[,-c(1,2)]

  drop <- caret::findLinearCombos(Tune_PRS)$remove
  drop <- names(Tune_PRS)[drop]

  Tune_PRS <- Tune_PRS %>% select(-all_of(drop))
  Validation_PRS <- Validation_PRS %>% select(-all_of(drop))

  if(family == "Gaussian"){
    Tune_PRS_tmp <- scale(Tune_PRS)
    xtx_1 <- apply(Tune_PRS_tmp,2,function(x){1/sum(x^2)})
    xty <- apply(Tune_PRS_tmp,2,function(x){sum(x*Tune_Y)})
    beta_hat <- xtx_1 * xty
    Predicted <- sweep(Tune_PRS_tmp, 2, beta_hat, "*")
    SSR <- colSums((Predicted - colMeans(Tune_Y))^2)
    SST <- sum((Tune_Y - colMeans(Tune_Y))^2)
    Order_Vec <- SSR/SST
  }else{
    Order_Vec <- apply(Tune_PRS,2,function(x){a <- data.frame(y = Tune_Y,x = x); b <- glm(y~x,data = a,family = binomial());return(auc(Tune_Y, unname(a$fitted.values)))})
  }

  prs_tune_order <- order(Order_Vec)
  Tune_PRS <- Tune_PRS[,prs_tune_order]
  Validation_PRS <- Validation_PRS[,prs_tune_order]

  tmp <- cor(Tune_PRS)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  #first remove prs (with smaller AUC) and have high corr with others
  Tune_PRS <- Tune_PRS[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
  Validation_PRS <- Validation_PRS[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]

  names_Tune_PRS <- colnames(Tune_PRS)
  names_Validation_PRS <- colnames(Validation_PRS)

  Tune_PRS <- data.frame(IID = Tune_IID,FID = Tune_FID,Tune_PRS)
  colnames(Tune_PRS) <- c("IID","FID",names_Tune_PRS)

  Validation_PRS <- data.frame(IID = Validation_IID,FID = Validation_FID,Validation_PRS)
  colnames(Validation_PRS) <- c("IID","FID",names_Validation_PRS)

  return(list(Cleaned_Tune_PRS = Tune_PRS,Cleaned_Validation_PRS = Validation_PRS))
}
