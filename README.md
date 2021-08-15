# missing_data_PCA
This project provides the R scripts (MissingData_PCA.r) and the simulated data sets used in the paper Yi & Latch 2021 "Nonrandom missing data can bias PCA inference of population genetic structure". Please cite this paper if reusing the R scripts or the simulated data sets. 

The R scripts included raw codes for generating coala models, simulating models and extract results into SNP matrices (function sim_raw_snp_matrix), incorporating missing data into SNP matrices (function add_missing_data), running PCA with the mean imputation approach in adegenet (function glPca), and plotting PCA sores. For a quick check of potential missing data effects, we recommend that researchers plot PCA scores with a color gradient showing per sample missingness. A cluster of high-missingness samples around the PCA origin, regardless of their genetic relationships, would indicate missing data effects and the samples around the PCA orign need to be interpreted with caution. 

