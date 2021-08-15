# missing_data_PCA
This project provides the R scripts (MissingData_PCA.r) and the simulated data sets used in the paper Yi & Latch 2021 "Nonrandom missing data can bias PCA inference of population genetic structure". Please cite this paper if reusing the R scripts or the simulated data sets. 

The R scripts included raw codes for generating coala models, simulating models and extract results into SNP matrices (function sim_raw_snp_matrix), incorporating missing data into SNP matrices (function add_missing_data), running PCA with the mean imputation approach in the package adegenet (function glPca; Jombart & Ahmed 2011), and plotting PCA sores using the package ggplot2 (Wickham 2016). 

For a quick check of potential missing data effects, we recommend that researchers plot PCA scores with a color gradient showing per sample missingness. A cluster of high-missingness samples around the PCA origin, regardless of their genetic relationships, would indicate missing data effects and the samples around the PCA orign need to be interpreted with caution. 


References 
Jombart T, Ahmed I. adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. 2011 Nov 1;27(21):3070-1. http://dx.doi.org/10.1093/bioinformatics/btr521
Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org
