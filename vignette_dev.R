
#install("~/Projects/DCEG/BB/CTSLEB_dev/CTSLEB")
library(CTSLEB)
library(data.table)

results.dir <- "test/"
data.dir <- "data/"
temp.dir <- "test/temp/"
system(paste0("mkdir -p ", temp.dir))
setwd('~/projects/CTSLEB_dev/')
plink19_exec <- "~/Apps/plink_v1.9/plink"
plink2_exec <- "~/Apps/plink2a/plink2"

# Step 1:

# Load summary statistics data for the reference (EUR) and the target (AFR)
# populations

sum_EUR <- fread(paste0(data.dir,"EUR_sumdata.txt"),header=T)
sum_AFR <- fread(paste0(data.dir,"AFR_sumdata.txt"),header=T)
head(sum_EUR)
head(sum_AFR)

# Step 1A: Alignment SNPs Across Populations

library(dplyr)

Eur_ref_plinkfile <- paste0(data.dir,"EUR_ref_chr22")
Afr_ref_plinkfile <- paste0(data.dir,"AFR_ref_chr22")
Afr_test_plinkfile <- paste0(data.dir,"AFR_test_chr22")
outprefix <- "chr22"
PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec)

dimCT(results_dir = results.dir,
      sum_target = sum_AFR,
      sum_ref = sum_EUR,
      ref_plink = Eur_ref_plinkfile,
      target_plink = Afr_ref_plinkfile,
      test_target_plink = Afr_test_plinkfile,
      out_prefix = outprefix,
      params_farm = PRS_farm)

pheno_file <- "data/y_tuning.txt"
y_tune <- fread(pheno_file)
PRS_tune(x = prs_mat,
         r2_vec = r2_vec,
         pthres = pthres,
         wc_base_vec = wc_base_vec,
         pheno_df = y_tune,
         tune_samples=10000,
         out = prs_out)
