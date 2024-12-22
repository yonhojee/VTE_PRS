# VTE_PRS
Code for Jee YH et al., Multi-ancestry polygenic risk scores for venous thromboembolism. Human Molecular Genetics 2024. https://doi.org/10.1093/hmg/ddae097. 
Data are publicly available at https://doi.org/10.7910/DVN/H3FR5U.

The LDPRED2 folder contained the following files: 
(1) a list of SNPs (info_snp_AFR) and their LDpred2-driven betas (beta_auto_AFR.Apr_03_2022) for AFR 
(2) a list of SNPs (info_snp_EUR) and their LDpred2-driven betas (beta_auto_EUR.Jun_10_2022) for EUR 
(3) a pseudo code that we used to generate the LDpred2_AFR, LDpred2_EUR. Please load info_snp and beta_auto files to run the code. 

The PRSCSX folder contained the following files: 
(1) 22 weight files for AFR (chr1 to chr22) [columns are: Chr, RSID, position, effect allele, Other allele, Weight] 
(2) 22 weight files for EUR (chr1 to chr22) [columns are: Chr, RSID, position, effect allele, Other allele, Weight] 
(3) a pseudo code that we used to generate the PRSCSx_AFR, PRSCSx_EUR using PLINK for your reference.
