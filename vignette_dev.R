#install("~/Projects/DCEG/BB/CTSLEB_dev/CTSLEB")
library(CTSLEB)
library(data.table)

results <- "test/"
data <- "data/"
temp <- "test/temp/"
system(paste0("mkdir -p ", temp))
setwd('~/projects/CTSLEB_dev/')
plink19_exec <- "~/Apps/plink_v1.9/plink"
plink2_exec <- "~/Apps/plink2a/plink2"

# STEP 1: Two-dimensional Clumping and Thresholding (CT)

# Load summary statistics data for the reference (EUR) and the target (AFR)
# populations

sum_EUR <- fread(paste0(data,"EUR_sumdata.txt"),header=T)
sum_AFR <- fread(paste0(data,"AFR_sumdata.txt"),header=T)
head(sum_EUR)
head(sum_AFR)

library(dplyr)

Eur_ref_plinkfile <- paste0(data,"EUR_ref_chr22")
Afr_ref_plinkfile <- paste0(data,"AFR_ref_chr22")
Afr_test_plinkfile <- paste0(data,"AFR_test_chr22")
outprefix <- "chr22"
PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec)

prs_mat <- dimCT(results_dir = results,
                 sum_target = sum_AFR,
                 sum_ref = sum_EUR,
                 ref_plink = Eur_ref_plinkfile,
                 target_plink = Afr_ref_plinkfile,
                 test_target_plink = Afr_test_plinkfile,
                 out_prefix = outprefix,
                 params_farm = PRS_farm)

pheno_tune_file <- "data/y_tuning.txt"
y_tune <- fread(pheno_tune_file)
prs_tune <- prs_mat[1:10000,]

n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test <- rep(0,n.total.prs)

for(p_ind in 1:n.total.prs){
  model <- lm(y_tune$V1~prs_tune[,(2+p_ind)])
  prs_r2_vec_test[p_ind] <- summary(model)$r.square
}

max_ind <- which.max(prs_r2_vec_test)
print(colnames(prs_tune)[max_ind+2])

# Step 2: Empirical-Bayes (EB) Estimation of Effect Sizes

prs_mat_eb <-CalculateEBEffectSize(bfile = Afr_test_plinkfile,
                              prs_tune = prs_tune,
                              plink_list = plink_list,
                              results_dir = results,
                              out_prefix = outprefix,
                              params_farm = PRS_farm)
