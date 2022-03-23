# CTSLEB
A multi-ancestry polygenic risk score approach


CTSLEB
=======
CT-SLEB is a method designed to generate multi-ancestry PRSs that incorporate existing large GWAS from EUR populations and smaller GWAS from non-EUR populations. The method has three key steps: 1. Clumping and Thresholding (CT) for selecting SNPs to be included in a PRS for the target population; 2. Empirical-Bayes (EB) method for estimating the coefficients of the SNPs; 3. Super-learning (SL) model to combine a series of PRSs generated under different SNP selection thresholds. The method requires three independent datasets: (1) GWAS summary statistics from training datasets across EUR and non-EUR populations; (2) a tuning dataset for the target population to find optimal model parameters; and (3) a validation dataset for the target population to report the final prediction performance. While this report assumes that individual-level data are available for model tuning and validation, summary-statistics-based methods could also be used in these steps. 

Usage
=======
[The 'CTSLEB' vignette](https://github.com/andrewhaoyu/TOP/blob/master/inst/TOP.pdf) will provide a good start point for using TOP package.


Development 
=======
This R package is developed by Haoyu Zhang and William Wheeler, and maintained by Haoyu Zhang <andrew.haoyu@gmail.com>.

Installation
=======
To install the development version of TOP, it's easiest to use the 'devtools' package.

install.packages("devtools")  
library(devtools)  
install_github("andrewhaoyu/TOP")

References
=======
Zhang, Haoyu, et al. "A mixed-model approach for powerful testing of genetic associations with cancer risk incorporating tumor characteristics." Biostatistics 22.4 (2021): 772-788.

