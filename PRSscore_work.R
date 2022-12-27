p_value_file_temp <- p_values
# k1 <- 5
# pthres[k1]
# idx <- which(unique_infor$P_ref<=pthres[k1])
temp.dir <- paste0(results.dir, "temp/")
for(k1 in 1:length(pthres)){
  #keep al the SNPs with P_EUR less than pthres[k1] in the analyses
  idx <- which(unique_infor$P_ref<=pthres[k1])
  p_value_file_temp$P[idx] = 0
  write.table(p_value_file_temp,file = paste0(temp.dir,"p_value_file"),col.names = F,row.names = F,quote=F)
  n_col <- ncol(scores)
  system(paste0(plink2_exec, " ",
                "--q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file ",
                "--score-col-nums 3-",n_col," ",
                "--score ",temp.dir,"score_file cols=+scoresums,-scoreavgs ",
                "--bfile ",data.dir,"AFR_test_chr22 ",
                "--out ",temp.dir,"prs_p_other_",k1))

}
