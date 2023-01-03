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
#system(paste0("mkdir -p ", temp))

plink19_exec <- "~/Apps/plink_v1.9/plink"
plink2_exec <- "~/Apps/plink2a/plink2"

Eur_ref_plinkfile <- paste0(data,"EUR_ref_chr22")
Afr_ref_plinkfile <- paste0(data,"AFR_ref_chr22")
Afr_test_plinkfile <- paste0(data,"AFR_test_chr22")
outprefix <- "chr22"
PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec)


load("step1_dev.RData")
#load("~/projects/CTSLEB_prod/step2_prod.RData")

# STEP2: Empirical-Bayes (EB) Estimation of Effect Sizes ####

best_snps <- colnames(prs_tune)[max_ind+2]
best_snp_set <- GetSNPSet(snp_ind = best_snps,
                          scores = scores,
                          clump_info = unique_infor)

unique_infor_post <- EBayesPostMean(clump_info = unique_infor,
                                    snp_set = best_snp_set)

# Get the posterior coefficient matrix

post_coef_mat <- cbind(unique_infor_post$BETA_EB_target,
                       unique_infor_post$BETA_EB_ref)
colnames(post_coef_mat) <- c("EB_target","EB_ref")

plink_list_eb <- PreparePlinkFileEBayes(snp_list = snp_list,
                                             clump_info = unique_infor,
                                             post_clump_info = unique_infor_post,
                                             post_beta = post_coef_mat,
                                             results_dir = results)
scores_eb <- plink_list_eb[[1]]
score_eb_file <- as.character(unlist(plink_list_eb["score_eb_file"]))
write.table(scores_eb,
            file = score_eb_file,
            row.names = F,
            col.names = F,
            quote=F)

p_values_eb <- plink_list_eb[[2]]

# Calculate the PRS based on the reference population Bayes coefficients and
# then combine all the prs

prs_mat_eb <- PRSscoreEbayes(bfile = Afr_test_plinkfile,
                             eb_plink_list = plink_list_eb,
                             results_dir = results,
                             out_prefix = outprefix,
                             params_farm = PRS_farm)

# We no have PRSs calculations based on EB coefficients for the reference and
# target population under all combinations of clumping r2-cutoff, window size,
# p-value thresholds. These PRS are the input for the super learning model which
# we will use to construct the PRS for the target population in STEP3
