library("devtools")
library(roxygen2)
library(data.table)
setwd('~/projects/CTSLEB_dev/')
#install("~/projects/CTSLEB_dev/CTSLEB")
library(CTSLEB)
library(dplyr)
dog_function()
dog_function(love=FALSE)

data <- "data/"
results <- "test/"
temp <- paste0(results,"temp/")
y_vad_file <- paste0(data,"y_validation.txt")
#system(paste0("mkdir -p ", temp))

plink19_exec <- "~/Apps/plink_v1.9/plink"
plink2_exec <- "~/Apps/plink2a/plink2"

load("step1_dev.RData")
load("step2_dev.RData")

#y_tune <- fread("data/y_tuning.txt")  # dataframe should already be present
y_vad <- fread(y_vad_file)

n.test <- 10000
prs_tune_sl <- as.data.frame(prs_mat_eb[1:n.test,-c(1:2),drop=F])
prs_vad <- as.data.frame(prs_mat_eb[(n.test+1):nrow(prs_mat_eb),-c(1:2),drop=F])

# drop all the prs columns with pairwise correlation more than 0.98

mtx <- cor(prs_tune_sl)

library(caret)

drop <- findCorrelation(mtx, cutoff=0.98)
drop <- names(prs_tune_sl)[drop]
prs_tune_sl_clean <- prs_tune_sl %>%
  select(-all_of(drop))
prs_vad_clean <- prs_vad %>%
  select(-all_of(drop))

#x <- prs_mat_eb

PRSTrainValSplit <- function(X,
                             n = 0.50) {
  mat_eb <- x
  n.test <- dim(mat_eb)[1]*n

  super_tune <- as.data.frame(mat_eb[1:n.test,-c(1:2),drop=F])
  super_validate <- as.data.frame(mat_eb[(n.test+1):nrow(mat_eb),-c(1:2),drop=F])

  # drop all the prs columns with pairwise correlation more than 0.98

  mtx <- cor(super_tune)

  drop <- findCorrelation(mtx, cutoff=0.98)
  drop <- names(super_tune)[drop]
  super_tune_clean <- super_tune %>%
    select(-all_of(drop))
  super_validate_clean <- super_validate %>%
    select(-all_of(drop))

  return_list <- list("tune" = super_tune_clean,
                      "validate" = super_validate_clean)
  return(return_list)
}

library(SuperLearner)
library(ranger)
library(glmnet)

train_val_list <- PRSTrainValSplit(X = prs_mat_eb)
prs_tune_sl_clean <- train_val_list[[1]]
prs_vad_clean <- train_val_list[[1]]


#choose the prediction algorithms

SL.libray <- c(
  "SL.glmnet",
  "SL.ridge",
  "SL.nnet"
  #"SL.bayesglm"
  #"SL.stepAIC"
  #"SL.xgboost"
  #"SL.randomForest"
  #"SL.ksvm",
  #"SL.bartMachine",
  #"SL.kernelKnn",
  #"SL.rpartPrune",
  #"SL.lm"
  #"SL.mean"
)

#train the super-learning model

sl <- SuperLearner(Y = y_tune$V1, X = prs_tune_sl_clean, family = gaussian(),
                  # For a real analysis we would use V = 10.
                  # V = 3,
                  SL.library = SL.libray)

#predict the outcome using the independent validation dataset

y_pred <- predict(sl, prs_vad_clean, onlySL = TRUE)

#evaluate the CT-SLEB prs performance on the validation

model <- lm(y_vad$V1~y_pred[[1]])
r2_ctsleb <- summary(model)$r.square
print(r2_ctsleb) ## [1] 0.000183132

TestModel <- function(x,y,

                      validate,
                      ml_library = c(
                        "SL.glmnet",
                        "SL.ridge",
                        "SL.nnet" )) {
  this_Y <- y
  this_X <- x

  sl <- SuperLearner(Y = this_Y, X = this_X, family = gaussian(),
                     SL.library = ml_libray)

  #predict the outcome using the independent validation dataset

  y_pred <- predict(sl, validate, onlySL = TRUE)

  #evaluate the CT-SLEB prs performance on the validation

  model <- lm(validate~y_pred[[1]])
}
