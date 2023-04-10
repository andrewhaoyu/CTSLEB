#' Create tune and validation data sets
#' @description. Rank the ebayes PRS based on R-square or AUC for the tuning
#' dataset and then remove highly correlated PRS from the tune and validate
#' mat_eb matrix.
#' @param x  tune subset of PRS matrix based on EB coefficients produced by
#' PRSscoreEBayes() or CalculateEBEffectSize().
#' @param y validate subset of PRS matrix based on EB coefficients produced by
#' PRSscoreEBayes() or CalculateEBEffectSize().
#' @param x_pheno vector containing tune phenotypes.
#' @param y_pheno vector containing validate phenotypes.
#' @param pheno_format 1 = continuous phenotype, 2 = binary phenotype. Default 1
#' @param params_farm List created by SetParamsFarm() function.
#' @return a list of dataframes for SuperLearner() function. [1] = tune_learn,
#' [2] = validate learn
#' @examples
#' sl_dataset <- Super_split(x = super_tune,
#' y = super_validate,
#' x_pheno = tune_pheno,
#' y_pheno = validate_pheno,
#' pheno_format = 2,
#' params_farm=PRS_farm)
#'
#' sl_tune_learn <- data.frame((sl_dataset[1]))
#' sl_val_learn <- data.frame(sl_dataset[2])
#'
#' @export

Super_split <- function(x,
                   y,
                   x_pheno,
                   y_pheno,
                   pheno_format = 1,
                   params_farm=as.null()){

  print("Executing Super_split() ... ")

  if (is.null(params_farm)) {
    print("no params_farm")
  } else {
    print("params_farm list will be used")
    plink19_exec <- as.character(unlist(params_farm["plink19_exec"]))
    plink2_exec <-  as.character(unlist(params_farm["plink2_exec"]))
    r2_vec <- as.numeric(unlist(params_farm["r2_vec"]))
    wc_base_vec <- as.integer(unlist(params_farm["wc_base_vec"]))
    memory <- as.integer(unlist(params_farm["mem"]))
    threads <- as.integer(unlist(params_farm["threads"]))
    pthres <- as.numeric(unlist(params_farm["pthres"]))
  }

  sl_tune <- x
  sl_val <- y
  sl_tune_pheno <- x_pheno
  sl_val_pheno <- y_pheno

  print("Executing correlation ... ")

  mtx <- cor(sl_tune)
  assign("mtx", mtx, envir = .GlobalEnv)

  qual_tune <- helper_qual_tune(x = sl_tune,
                            x_pheno = sl_tune_pheno,
                            pthres = pthres,
                            r2_vec= r2_vec,
                            wc_base_vec = wc_base_vec,
                            pheno_format = pheno_format)

  ix_keep <- helper_qual_order(x = mtx, y = qual_tune)
  prs_tune_learn <- data.frame(sl_tune[,ix_keep])
  prs_val_learn <- data.frame(sl_val[,ix_keep])

  print("tune object created ... ")
  print("validate object created ... ")

  return_list <- list(prs_tune_learn,prs_val_learn)

  return(return_list)
}
