results.dir <- "test/"
data.dir = "data/"
temp.dir <- "test/temp/"
system("mkdir test/temp/")

plink19_exec <- "/Users/dayneokuhara/Apps/plink_apps/plink1.9/plink"
plink2_exec <- "/Users/dayneokuhara/Apps/plink_apps/plink2"
#install_github("CTSLEB")

library(CTSLEB)
library(data.table)

# Step 1:

# Load summary statistics data for the reference (EUR) and the target (AFR)
# populations

sum_EUR <- fread(paste0(data.dir,"EUR_sumdata.txt"),header=T)
sum_AFR <- fread(paste0(data.dir,"AFR_sumdata.txt"),header=T)
head(sum_EUR)
head(sum_AFR)

# Step 1A: Alignment SNPs Across Populations

library(dplyr)

# Align ref GWAS summary stat to the target population and then split based on
# p-values

AlignSum(sum_target = sum_AFR,
         sum_ref = sum_EUR,
         results_dir = results.dir,
         SplitSum = FALSE)
head(sum_com)
sum_com_split <- SplitSum(sum_com)

WriteSplitTables(sum_com_split,results_dir = results.dir)

# Alternatively, setting SplitSum = TRUE will automatically perform
# the SplitSum() and WriteSplitTables() functions and create output variable
# and file names

AlignSum(sum_target = sum_AFR,
         sum_ref = sum_EUR,
         results_dir = results.dir,
         SplitSum = TRUE)

# Step 1B: Clumping with Plink1.9

Eur_plinkfile <- paste0(data.dir,"EUR_ref_chr22")
Afr_plinkfile <- paste0(data.dir,"AFR_ref_chr22")
outprefix <- "chr22"

# default settings for RunClump()
# r2_vec <- c(0.01,0.05,0.1,0.2,0.5,0.8)
# wc_base_vec <- c(50,100)
# mem = 8000
# threads = 2

RunClump(plink19_exec,
         ref_plink = Eur_plinkfile,
         target_plink = Afr_plinkfile,
         out_prefix = outprefix,
         results_dir = results.dir
         )

# Alternatively, both the AlignSum() and RunClump() can be performed using a
# clump_farm object which is simply a list of options passed to
# the AlignSum() and RunClump() 

PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec)
dimCT(plink19_exec,
      results_dir = results.dir,
      sum_target = sum_AFR,
      sum_ref = sum_EUR,
      ref_plink = Eur_plinkfile,
      target_plink = Afr_plinkfile,
      out_prefix = outprefix,
      params_farm = PRS_farm)

# Prepare plink2 files for PRS calculations

# plink2 requires three files to compute PRS under different thresholds
# The first file is a score_file, which contains the SNP coefficients
# The second file is a p_value_file, which contains the SNP p-values
# The third file is a q_range file, which contains the p_value thresholds
# PreparePlinkile() will prepare these three files and if return_list is TRUE
# create a global variable 'plink_files' which contains these four data frames
# in a list.

#pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)

PreparePlinkFile(snp_list = snp_list,
                 sum_com = sum_com,
                 results_dir = results.dir,
                 return_list = TRUE)
class(plink_files)
names(plink_files)
scores <- plink_files[[1]]
scores_file <- paste0(temp.dir,"score_file")
write.table(scores,
            file = scores_file ,
            row.names = F,
            col.names = F,
            quote=F)

p_values <- plink_files[[2]]
unique_infor <- plink_files[[3]]

q_range <- plink_files[[4]]
head(q_range)
q_range_file <- paste0(temp.dir,"q_range_file")
write.table(q_range,
            file = q_range_file,
            row.names = F,
            col.names = F,
            quote=F)

# Alternative: Use the return_list=FALSE default to automatically created the
# the global variables q_range, score_file, p_value_file and unique_infor

rm(q_range)
rm(scores)
rm(scores_file)
rm(q_range_file)
rm(unique_infor)
rm(plink_files)

PreparePlinkFile(snp_list = snp_list,
                 sum_com = sum_com,
                 results_dir = results.dir,
                 return_list = FALSE)

# Generate plink2 PRS
# using default pthers values

PRSscore(plink2_exec = plink2_exec,
         bfile = Afr_plinkfile,
         scores = scores,
         p_values = p_values,
         out_prefix = outprefix,
         threads = 4,
         memory = 8000,
         results_dir = results.dir
)

# Combine plink2 PRS results

prs_list <- list()
temp <- 1
#take the column name of different clumping parameters
names <- colnames(scores)[3:ncol(scores)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    #the --score file cols=+scoresums,-scoreavgs command in plink2 computes PRS as G*beta
    #If you compute PRS by chromosome, you need to sum the PRS scores for all chromosomes.
    #load PRS for SNPs with p_target<p_thres[k2]|p_eur<p_thres[k1]
    prs_temp <- fread(paste0(temp.dir,"prs_p_other_",k1,".p_tar_",k2,".sscore"))
    # times (2*number of SNPs)
    prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
    colnames(prs_list[[temp]]) = paste0(names,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp <- temp + 1
  }
}

prs_mat <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))

# Alternatively, both the PRS generation with plink2 and combining the results
# into a single matrix can be accomplished with the PRSscore() function and
# PRS_farm list

prs_mat <- PRSscore(plink2_exec = plink2_exec,
                    bfile = AFR_plinkfile,
                    q_range_file = q_range_file,
                    scores = scores,
                    scores_file = scores_file,
                    p_values = p_values,
                    pthres = pthres,
                    threads = 4,
                    memory = 8000)

# Selecting the PRS best snps

# take the first 10,000 subjects for tuning purpose
prs_tune <- prs_mat[1:10000,]
# find the best R-square among the all the PRSs to find candidate set
# we use this candidate for estimating covariance matrix for the prior distribution #create prediction r2 vector to store r2 for different prs
n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test <- rep(0,n.total.prs)
#load the phenotype data for the tuning set
y_tune <- fread("data/y_tuning.txt")
for(p_ind in 1:n.total.prs){
  #the first two columns of prs_tun are family id and individual id #prs starts from the third column
  model = lm(y_tune$V1~prs_tune[,(2+p_ind)])
  prs_r2_vec_test[p_ind] = summary(model)$r.square
}
max_ind <- which.max(prs_r2_vec_test)
#+2 is due to the first two columns are family id and individual id print(colnames(prs_tun)[max_ind+2])
print(colnames(prs_tun)[max_ind+2])

# Alternative

pheno_file <- "data/y_tuning.txt"
y_tune <- fread(pheno_file)
PRS_tune(x = prs_mat,
         r2_vec = r2_vec,
         pthres = pthres,
         wc_base_vec = wc_base_vec,
         pheno_df = y_tune,
         tune_samples=10000,
         out = prs_out)
