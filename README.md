# CTSLEB
A multi-ancestry polygenic risk score approach


CTSLEB
=======
CT-SLEB is a method designed to generate multi-ancestry PRSs that incorporate existing large GWAS from EUR populations and smaller GWAS from non-EUR populations. The method has three key steps: 1. Clumping and Thresholding (CT) for selecting SNPs to be included in a PRS for the target population; 2. Empirical-Bayes (EB) method for estimating the coefficients of the SNPs; 3. Super-learning (SL) model to combine a series of PRSs generated under different SNP selection thresholds. The method requires three independent datasets: (1) GWAS summary statistics from training datasets across EUR and non-EUR populations; (2) a tuning dataset for the target population to find optimal model parameters; and (3) a validation dataset for the target population to report the final prediction performance. While this report assumes that individual-level data are available for model tuning and validation, summary-statistics-based methods could also be used in these steps. 

Usage
=======
[The 'CTSLEB' vignette](https://andrewhaoyu.github.io/CTSLEB/) will provide a good start point for using CTSLEB package.

Development 
=======
This R package is developed by Haoyu Zhang, and maintained by Haoyu Zhang <andrew.haoyu@gmail.com>.

Installation
=======
To install CTSLEB, it's easiest to use the 'devtools' package.

install.packages("devtools")  
library(devtools)  
install_github("andrewhaoyu/CTLEB")

Data
=======
The example data for the tutorial can be downloaded through this [link](https://drive.google.com/file/d/1wswLKQmgYgkkog_vaDaVlLEmgoQS_xLG/view?usp=sharing).

For data analyses, the pipeline require reference sample from different populations for the clumping step. We use data from 1000 Genomes Project (Phase 3) in PLINK format. Other reference data is also fine as long as it's in PLINK format. The 1000 Genomes Project data in PLINK format can be download through the following links:
[African (AFR) populations]()
[American (AMR) populations]()
[European (EUR) populations]()
[East Asian (EAS) populations]()
[South Asian (SAS) populations]()

Other software
=======
The analyses also need [PLINK 1.9](https://www.cog-genomics.org/plink/) for clumping purpose and [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) for calculating PRSs. Guidance can be found in the vignette.

Reference
=======
Zhang, Haoyu, et al. "Novel Methods for Multi-ancestry Polygenic Prediction and their Evaluations in 3.7 Million Individuals of Diverse Ancestry" Biorxiv (2022).

