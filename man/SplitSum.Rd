% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SplitSum.R
\name{SplitSum}
\alias{SplitSum}
\title{Inputs an AlignSum() object and splits SNPs into two groups}
\usage{
SplitSum(x, results_dir, write_tables = TRUE)
}
\arguments{
\item{x}{sum_com Data.frame object produced by AlignSum().}

\item{results_dir}{Folder or directory to export files. Mandatory when
'write_tables = TRUE'. Default file names are sum_ref.txt and sum_target.txt.}

\item{write_tables}{Writes two files to the 'results_dir' used by the plink1.9
'--clump' parameter. Default is TRUE.}
}
\value{
A list with the two data.frame objects.SNPs with p-values in the
ref population smaller than the target population in the [[1]] slot.
SNPs with p-values in the target population smaller than the ref population
in the [[2]] slot. If "write_tables = TRUE" the tables will be written as files
in the results_dir folder AND the global variables ref_split_file and
target_split_file will be created for the file locations.
}
\description{
Split the SNPs into two groups: 1. SNPs with p-values in the ref population
smaller than the target population. 2. SNPs with p-values in the target
population smaller than the ref population.
}
\examples{
data <- "data/"
results <- "test/"
sum_EUR <- fread(paste0(data, "EUR_sumdata.txt"), header=T)
sum_AFR <- fread(paste0(data, "AFR_sumdata.txt"), header=T)
sum_com <- AlignSum(sum_tar = sum_AFR, sum_other = sum_EUR)
split_list <- SplitSum(sum_com = sum_com,
                       results_dir = results,
                       write_tables = TRUE)
}
