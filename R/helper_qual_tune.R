#' Helper function for Super_split. Create qual_tune vector consisting
#' of either r2 square for linear regression (pheno_format = 1)
#' or AUC for binary (pheno_format = 2)
#'
#' @param x tune matrix PRS values
#' @param x_pheno vector of tune phenotype values
#' @param pthres vector of tune phenotype values
#' @param r2_vec vector of tune phenotype values
#' @param wc_base_vec vector of tune phenotype values
#' @param r2_vec vector of tune phenotype values
#' @return vector

helper_qual_tune <- function(x,
                             x_pheno,
                             pthres,
                             y,
                             y_pheno,
                             r2_vec,
                             wc_base_vec,
                             pheno_format = 1) {

  if (pheno_format == 1) {
    print("Executing helper_qual_tune() for continuous phenotype... ")
    sl_tune <- x
    sl_tune_pheno <- x_pheno
    n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
    qual_tune <- rep(0,n.total.prs)
    for(p_ind in 1:n.total.prs){
      model <- lm(sl_tune_pheno~sl_tune[,(2+p_ind)])
      #print(paste0(colnames(sl_tune)[2+p_ind],
      #             ": qual = ",
      #             summary(model)$r.square))
      qual_tune[p_ind] <- summary(model)$r.square
    }
    print(paste0("max tune quality: ", max(qual_tune)))
    print(paste0("min tune quality: ", min(qual_tune)))
    return(qual_tune)
  } else if (pheno_format == 2) {
    print("Executing helper_qual_tune() for binary phenotype... ")
    sl_tune <- x
    sl_tune_pheno <- x_pheno
    sl_val <- x
    sl_val_pheno <- x_pheno
    auc_tune <- rep(0,n.total.prs)
    for(p_ind in 1:n.total.prs){
      tune <- data.frame(V1=sl_tune_pheno,V2=sl_tune[,(2+p_ind)])
      val <- data.frame(V1=sl_val_pheno,V2=sl_val[,(2+p_ind)])
      model <- glm(V1~V2, data = tune, family="binomial")
      prediction <- predict(model, newdata=val, type="response")
      auc_tune[p_ind] <- auc(val$V1, prediction, quiet=TRUE)[1]
    }
    print(paste0("max tune quality: ", max(auc_tune)))
    print(paste0("min tune quality: ", min(auc_tune)))
    return(auc_tune)
  } else {
    stop('"please select either 1 (continuous) or 2 (binary) for "pheno_format"')
  }
}

