CTSLEB
=======
CT-SLEB is a method designed to generate multi-ancestry PRSs that incorporate existing large GWAS from EUR populations and smaller GWAS from non-EUR populations. The method has three key steps: 1. Clumping and Thresholding for selecting SNPs to be included in a PRS for the target population; 2. Empirical-Bayes method for estimating the coefficients of the SNPs; 3. Super-learning model to combine a series of PRSs generated under different SNP selection thresholds. The method requires three independent datasets: (1) GWAS summary statistics from training datasets across EUR and non-EUR populations; (2) a tuning dataset for the target population to find optimal model parameters; and (3) a validation dataset for the target population to report the final prediction performance. 

Usage
=======
[The 'CTSLEB' vignette](https://andrewhaoyu.github.io/CTSLEB/) will provide a good start point for using CTSLEB package.

Installation
=======
To install CTSLEB, it's easiest to use the 'devtools' package.

install.packages("devtools")  
library(devtools)  
install_github("andrewhaoyu/CTSLEB")

Data
=======
The example dataset for the vignette can be downloaded through this [link](https://drive.google.com/file/d/1wswLKQmgYgkkog_vaDaVlLEmgoQS_xLG/view?usp=sharing).

For data analyses, the pipelines require reference samples from different populations for the clumping step. We use data from 1000 Genomes Project (Phase 3) in PLINK format. Other reference data can also be used as long as it's in PLINK format. The 1000 Genomes Project data in PLINK format can be downloaded through the following links:

[African (AFR) populations](https://drive.google.com/file/d/1pwQuM3pb8qSNyNyEsCeoy_jqz5CLIjg_/view?usp=sharing)

[American (AMR) populations](https://drive.google.com/file/d/1xDCPvj-JPTbWOS9fm3jovLrFoXJmpOb2/view?usp=sharing)

[European (EUR) populations](https://drive.google.com/file/d/1cWjUfDfar-shXbyLvbDfW-BqCZ-yP6h4/view?usp=sharing)

[East Asian (EAS) populations](https://drive.google.com/file/d/1xrkzc06RG6KcjYctW9VDG9VnQj2Odczh/view?usp=sharing)

[South Asian (SAS) populations](https://drive.google.com/file/d/1n2yR2ZoMHCS_UHr2wMzACXF9wDWCVu6r/view?usp=sharing)

Data sharing
=======
Simulation can be challenging and time-consuming to generate large-scale genotype data with realistic LD for diverse ancestries. We have generated independent 600,000 subjects for five ancestries: African, American, European, East Asian, and South Asian using the LD estimated from the 1000 Genomes Project (Phase 3). Each of the five ancestries contains 120,000 subjects. We release this data through Harvard DataVerse to help researchers develop and test methods in multi-ancestry GWAS settings. The data can be downloaded through the following [link](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP). All the simulated genotype data, phenotype data, and summary statistics used in the manuscript are available. The ReadMe file in the directory gives a detailed explanation of the data. Due to the file size limit on Harvard DataVerse, we restrict the shared data on SNPs on HapMap3 (HM3) or Multi-Ethnic Genotyping Arrays (MEGA) chip arrays, or both. If you need the data for all ~19 million SNPs, please contact Haoyu Zhang (andrew.haoyu@gmail.com).

Other software
=======
The analyses also need [PLINK 1.9](https://www.cog-genomics.org/plink/) for clumping purpose and [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) for calculating PRSs. Guidance can be found in the vignette.

Support 
=======
Please direct any problems or questions to Haoyu Zhang <andrew.haoyu@gmail.com>.

Reference
=======
Zhang, Haoyu, et al. "Novel Methods for Multi-ancestry Polygenic Prediction and their Evaluations in 3.7 Million Individuals of Diverse Ancestry" Biorxiv (2022).

