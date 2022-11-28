
system("mkdir test/temp.dir")
temp.dir = "test/temp.dir/"
soft.dir = "test/"
data.dir = "data/"
plink19_exec <- "/Users/dayneokuhara/Apps/plink_apps/plink1.9/plink"

#install_github("CTSLEB")

library(CTSLEB)
library(data.table)

# Load summary statistics data for the reference (EUR) and the target (AFR)
# populations

sum_EUR <- fread(paste0(data.dir,"EUR_sumdata.txt"),header=T)
sum_AFR <- fread(paste0(data.dir,"AFR_sumdata.txt"),header=T)
head(sum_EUR)
head(sum_AFR)

# Step 1: Alignment SNPs Across Populations

library(dplyr)

# Align ref GWAS summary stat to the target population and then split based on
# p-values

AlignSum(sum_target = sum_AFR,
         sum_ref = sum_EUR,
         SplitSum = FALSE)
head(sum_com)
sum_com_split <- SplitSum(sum_com)

ref_split_file <- paste0(temp.dir,"sum_EUR.txt")
target_split_file <- paste0(temp.dir,"sum_AFR.txt")

WriteSplitTables(sum_com_split,
                 ref_split_file = ref_split_file,
                 target_split_file = target_split_file)

# Alternatively, setting SplitSum = TRUE will automatically perform
# the SplitSum() and WriteSplitTables() functions with default output variable
# and file names

AlignSum(sum_target = sum_AFR,
         sum_ref = sum_EUR,
         ref_split_file = ref_split_file,
         target_split_file = target_split_file,
         SplitSum = TRUE)

# Clumping with Plink1.9
# Plink1.9 is used to ld clump the reference and target snps based on a range of
# r2 and base window sizes. The clumping_window_size = base_window_size/clumping_r_square
# so lower r2 thresholds have a larger clumping window size.

Eur_plinkfile <- paste0(data.dir,"EUR_ref_chr22")
Afr_plinkfile <- paste0(data.dir,"AFR_ref_chr22")

r2_vec <- c(0.01,0.05,0.1,0.2,0.5,0.8)
wc_base_vec <- c(50,100)

snp_list <- list()
temp <- 1
for(r_ind in 1:length(r2_vec)){
  wc_vec = wc_base_vec/r2_vec[r_ind]
  for(w_ind in 1:length(wc_vec)){
    pthr <- 1
    r2thr <- r2_vec[r_ind]
    kbpthr <- wc_vec[w_ind]
    ref_outfile <- paste0(temp.dir,"ref_CT_rind_",r_ind,"_wcind_",w_ind)
    print(ref_outfile)
    target_outfile <- paste0(temp.dir,"target_CT_rind_",r_ind,"_wcind_",w_ind)
    Plink19Clump(plink19_exec,
                 bfile = Eur_plinkfile,
                 clump = ref_split_file,
                 clump_p1 = pthr,
                 clump_r2 = r2thr,
                 clump_kb = kbpthr,
                 out = ref_outfile)
    Plink19Clump(plink19_exec,
                 bfile = Afr_plinkfile,
                 clump = target_split_file,
                 clump_p1 = pthr,
                 clump_r2 = r2thr,
                 clump_kb = kbpthr,
                 out = target_outfile)

    LD_ref <- fread(paste0(ref_outfile,".clumped"))[,3,drop=F]
    LD_target<- fread(paste0(target_outfile,".clumped"))[,3,drop=F]
    LD  <- rbind(LD_ref,LD_target)
    snp_list[[temp]] <- LD
    names(snp_list[[temp]]) <- paste0("clump_r2_",r2thr,"_ws_",kbpthr)
    temp <- temp + 1
  }
}

# Alternatively, both the AlignSum() and Plink19Clump() can be performed using a
# clump_farm object which is simply a list of options which can be passed to
# the AlignSum() and Plink19Clump()

ref_out <- "test/temp/ref"
target_out <- "test/temp/target"
clump_farm <- SetClumpFarm()
dimCT(plink19_exec,
      sum_target = sum_AFR,
      sum_ref = sum_EUR,
      ref_split_file = ref_split_file,
      target_split_file = target_split_file,
      ref_plink = Eur_plinkfile,
      target_plink = Afr_plinkfile,
      ref_clump_out = ref_out,
      target_clump_out = target_out,
      clump_farm=clump_farm)

# Thresholding with Plink2a_dev

# plink2 requires three files to compute PRS under different thresholds
# The first file is a score_file, which contains the SNP coefficients
# The second file is a p_value_file, which contains the SNP p-values
# The third file is a q_range file, which contains the p_value thresholds
# PreparePlinkile() will prepare these three files and if return_list is TRUE
# create a global variable 'plink_files' which contains these four data frames
# in a list.

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
PreparePlinkFile(snp_list,
                 sum_com,
                 pthres = pthres,
                 output = temp.dir,
                 return_list = TRUE)
score_file <- plink_files[[1]]
head(score_file)
p_value_file <- plink_files[[2]]
head(p_value_file)
unique_infor <- plink_files[[3]]
head(unique_infor)
q_range <- plink_files[[4]]
head(q_range)

write.table(score_file,file = paste0(temp.dir,"score_file"),row.names = F,col.names = F,quote=F)
write.table(q_range,file = paste0(temp.dir,"q_range_file"),row.names = F,col.names = F,quote=F)

# Alternative: Use the return_list=FALSE default to automatically created the
# the global variables q_range, score_file, p_value_file and unique_infor

rm(q_range)
rm(score_file)
rm(p_value_file)
rm(unique_infor)
rm(plink_files)

PreparePlinkFile(snp_list,
                 sum_com,
                 pthres = pthres,
                 output = temp.dir,
                 return_list = FALSE)
