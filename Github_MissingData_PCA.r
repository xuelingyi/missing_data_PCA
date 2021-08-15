## This script includes the R codes used in the paper Yi & Latch 2021 "Nonrandom missing data can bias PCA inference of populaiton genetic structure". 

library(coala)
activate_ms()

library(adegenet)
library(stringr)
library(ggplot2)
library(ggpubr)

####### define evolutionary models in coala ###########
nloci = 5000 
nindv = 25

models=c("p3", "p3_mig50", "island_mig50_mig5", "cline_mig50")
n_pops = c(3, 3, 3, 3, 3)

p3 = coal_model(sample_size=c(rep(nindv, 3)), loci_number=nloci, loci_length=1, ploidy=2) +
  feat_mutation(rate = 1, fixed_number = T) +
  feat_pop_merge(0.9, 3,2) + feat_pop_merge(0.9, 2,1) +
  sumstat_seg_sites()
##  divergence of all populations at 0.9 time units backwards
## with the default "IFS" mutation model, but probably does not matter because the length is fixed to 1

p3_mig50 = p3 + feat_migration(rate = 50, symmetric = T)
## the feat_migration_rate is 4Nm, and N=nindv=25, meaning that the migration rate m=0.5 in this model

cline_mig50 = p3 + 
  feat_migration(rate = 50, 1,2) + feat_migration(rate = 50, 2,1) +
  feat_migration(rate = 50, 2,3) + feat_migration(rate = 50, 3,2) 

island_mig50_mig5 = p3 + 
  feat_migration(rate = 50, 1,3) + feat_migration(rate = 50, 3,1) +
  feat_migration(rate = 5, 1,2) + feat_migration(rate = 5, 2,1) +
  feat_migration(rate = 5, 3,2) + feat_migration(rate = 5, 2,3)
  
####### simulate coala models and extract population genetic data sets as SNP matrices #########
sim_raw_snp_matrix = function(coala_model, pop_size, individual_per_population, loci_number, ...){
    
	myrownames = paste(rep(paste0("pop", 1:pop_size), each=individual_per_population), 
                       rep(paste0("indv", 1:individual_per_population),  pop_size), sep=".")
	
	sim = simulate(get(coala_model), seed=sample(1:500,1))
    data = unlist(sim$seg_sites, recursive = F)
    names(data) = paste0(rep(c("snps", "position", "trio_locus"), loci_number), 
	                     rep(1:loci_number, each=3))
    snps = as.data.frame(data[paste0("snps", 1:loci_number)])
    SNP = snps[1,] + snps[2,]
    for (i in 2:(pop_size*individual_per_population)) {
      snp_indv = snps[2*i-1, ] + snps[2*i, ]
      SNP = rbind(SNP, snp_indv)
    }
    row.names(SNP) = myrownames
	
	return(SNP)
}

for (m in 1:length(models)) {
  
  npop = n_pops[m]
  SNP = sim_raw_snp_matrix(models[m], npop, nindv, nloci)  
  save(npop, SNP, file=paste0("./simulation/", models[m], "_SNP.RData"))
}

# replicate simulates
for (m in 1:length(models)) {
  
  for (rep in 1:4) {
  npop = n_pops[m]
  SNP = sim_raw_snp_matrix(models[m], npop, nindv, nloci)  
  save(npop, SNP, file=paste0("./simulation/", models[m], "_replicate", rep, "_SNP.RData"))
  }
}


###### introduce missing data into the raw SNP matrices ########
add_missing_data = function(SNP_matrix, #the simulated raw SNP matrix
                            indv_per_pop, loci_number, pop_size,
                            total_miss, #average percent of missing data per population
                            random=T, #generate data frame "rand"
                            bias_indv=T, #generate data frame "biasINDV"
                            bias_pop=T, #generate data frame "biasPOP"
                            miss_bias=0.8, #80% missing data are biased
                            miss_bias_indv=0.2, #biased missing data in 20% individuals
                            miss_bias_popID=2, #biased missing data in pop2
                            ...) { 
  
  # amount of missing data per population (n) and total (N)
  n = indv_per_pop*loci_number*total_miss
  N = n*pop_size
  
  rand = NULL 
  biasINDV = NULL 
  biasPOP= NULL 
  biasINDVPOP = NULL
  
  ## random missing data
  if (random == T) {
    rand = SNP_matrix
    
    rand = unlist(rand)
    rand[sample(length(rand), N)] = NA
    rand = as.data.frame(matrix(rand, ncol=loci_number))
    row.names(rand) = rownames(SNP_matrix)
	
	rand <<- rand
  }
  
  ## individual-biased missing data
  if (bias_indv == T){
    biasINDV = SNP_matrix
    
    for (p in 1:pop_size) {
    #default 80% (miss_bias) missing data condensed in 20% (miss_bias_indv) individuals per population
      indID = sample(indv_per_pop, indv_per_pop*miss_bias_indv)
      
      pop = biasINDV[((p-1)*indv_per_pop+1):(p*indv_per_pop),]
      
      biasindv = unlist(pop[indID, ])
      biasindv[sample(length(biasindv), n*miss_bias)] = NA
      biasindv = as.data.frame(matrix(biasindv, ncol=loci_number))
      
      otherindv = unlist(pop[-indID, ])
      otherindv[sample(length(otherindv), n*(1-miss_bias))] = NA
      otherindv = as.data.frame(matrix(otherindv, ncol=loci_number))
      
	  biasINDV[((p-1)*indv_per_pop+1):(p*indv_per_pop),] = rbind(biasindv, otherindv)
	  
	  if (p==1) {newrownames = c(rownames(pop[indID, ]), rownames(pop[-indID, ]))} else {newrownames = c(newrownames, rownames(pop[indID, ]), rownames(pop[-indID, ]))}
	   
    }
   rownames(biasINDV) = newrownames
   biasINDV <<- biasINDV
  }
  
  ## population-biased missing data
  if (bias_pop == T){
    biasPOP = SNP_matrix
    
    for (p in 1:pop_size) {
      # default 80% (miss_bias) missing data in pop2 (miss_bias_popID) individuals (random within pop2)
      pop = biasPOP[((p-1)*indv_per_pop+1):(p*indv_per_pop),]
      
      pop = unlist(pop)
      pop[sample(length(pop), 
                 ifelse(p == miss_bias_popID, N*miss_bias, N*(1-miss_bias)/(pop_size-1)))] = NA
      pop = as.data.frame(matrix(pop, ncol=loci_number))
      biasPOP[((p-1)*indv_per_pop+1):(p*indv_per_pop),] = pop
    }  
    biasPOP <<- biasPOP
    }

}

for (m in 1:length(models)) {
  print(models[m])
  load(paste0("./simulation/", models[m], "_SNP.RData"))
  # loaded data: npop, SNP
  
  for (miss in c(0.01, 0.1, 0.2)) { 
    add_missing_data(SNP, nindv, nloci, npop, total_miss=miss, 
                     random=T, bias_indv=T, bias_pop=T, 
                     miss_bias=0.8, miss_bias_indv=0.2, miss_bias_popID=2)
	
	save(npop, miss, rand, 
         file=paste0("./MISSdata/", models[m], "_SNP_rand_miss", miss, ".RData"))
    save(npop, miss, biasINDV, 
         file=paste0("./MISSdata/", models[m], "_SNP_biasINDV_miss", miss, ".RData"))
    save(npop, miss, biasPOP, 
         file=paste0("./MISSdata/", models[m], "_SNP_biasPOP_miss", miss, ".RData"))
    
	rm(rand, biasINDV, biasPOP)
  }
  rm(npop, SNP)
}

## in the cline modle: change the biased population from pop2 (the admixed population) into pop3
model="cline_mig50"
load(paste0("./simulation/", model, "_SNP.RData"))
# loaded data: npop, SNP
for (miss in c(0.01, 0.1, 0.2)) { 
  add_missing_data(SNP, nindv, nloci, npop, total_miss=miss, 
                   random=F, bias_indv=F, bias_pop=T, 
				   miss_bias=0.8, miss_bias_indv=0.2, 
				   miss_bias_popID=3)
				   ## Change the biased population into pop3
	
  save(npop, miss, biasPOP, 
       file=paste0("./MISSdata/", model, "_SNP_biasPOP_pop3_miss", miss, ".RData"))
  rm(biasPOP)
}
rm (npop, SNP, model, miss)

# can double check total percent of missing data: should be equal to miss
##  length(unlist(rand)[is.na(unlist(rand))])/(nloci*npop*nindv)
# can double check relative percent of missing data in pop2: should be 0.8
##  length(unlist(biasPOP[(nindv+1):(nindv*2),])[is.na(unlist(biasPOP[(nindv+1):(nindv*2),]))])/length(unlist(biasPOP)[is.na(unlist(biasPOP))])
# can double check relative percent of missing data in biased individuals: should be 0.8
##  length(unlist(biasINDV[(nindv+1):(nindv*1.2),])[is.na(unlist(biasINDV[(nindv+1):(nindv*1.2),]))])/length(unlist(biasINDV[(nindv+1):(nindv*2),])[is.na(unlist(biasINDV[(nindv+1):(nindv*2),]))])


######## run PCA in adegenet (default mean imputation on missing data) #############
run_glPCA = function(SNP_data, loci_number,...) {
  
  gl = new("genlight", SNP_data, ploidy=2)
  pca = glPca(gl, nf=loci_number) # nloci > nindv, keep all axes
  scores = as.data.frame(pca$scores)
  scores$indv = rownames(scores)
  scores$pop = gsub(".indv.*", "", scores$indv)
  
  indvmiss = NULL
  for (i in 1:nrow(SNP_data)) {
    x = SNP_data[i, ]
    y = length(x[x == "NA"]) / loci_number
    indvmiss = c(indvmiss, y)
  }
  scores$indvmiss = round(indvmiss, 2)
  
  snpmiss = NULL
  for (i in 1:loci_number) {
    x = SNP_data[, i]
    y = length(x[x == "NA"]) / nrow(SNP_data)
    snpmiss = c(snpmiss, y)
  }
  
  pca <<- pca
  scores <<- scores
  snpmiss <<- snpmiss

}

## raw SNP matrices in ./simulation/
raw_data = list.files("./simulation")
for (m in 1:length(raw_data)) {
  
  load(paste0("./simulation/", raw_data[m]))
  # loaded data: npop, SNP 
  run_glPCA(SNP, nloci)
  save(pca, scores, snpmiss, file=paste0("./PCAdata/", raw_data[m], "_PCA.RData"))
  
  rm(SNP, npop, pca, scores, snpmiss)
  rm(miss, npop, pca, scores, snpmiss)
}
rm(raw_data)

## SNP matrices with missing data in ./MISSdata/
miss_data = list.files("./MISSdata/")
for (m in 1:length(miss_data)) {
  
  load(paste0("./MISSdata/", miss_data[m]))
  # loaded data: npop, gsub(".*SNP_|_miss.*", "", miss_data[m]) --> rand or biasINDV or biasPOP, miss
  run_glPCA(get(gsub(".*SNP_|_miss.*", "", miss_data[m])), nloci)
  save(pca, scores, snpmiss, file=paste0("./PCAdata/", miss_data[m], "_PCA.RData"))
  
  rm(list=gsub(".*SNP_|_miss.*", "", miss_data[m]))
  rm(miss, npop, pca, scores, snpmiss)
}
rm(miss_data)

######### plot PCA results with color coded missingness per individual #########
plot_PCA_scores = function(data_name, scores_file, pca_file, my_shape, print=T, plot_title=NULL, shape_guide=T, ...){

  myplot = ggplot(scores_file[order(scores_file$indvmiss, decreasing=F),], 
                  aes(x=PC1, y=PC2, shape=pop, fill=indvmiss)) +       #### NOTE: fill color as a gradient based on individual missing values in the input data set; recommended for a quick check of potential missing data effects 
    geom_point(size=1.8, stroke=0.02, alpha=0.9) + 
    scale_shape_manual(values=my_shape) +
    guides(shape=ifelse(shape_guide==T, "legend", F))+
    geom_hline(yintercept=0, color="grey60", size=0.1) + 
    geom_vline(xintercept=0, color="grey60", size=0.1) +
    theme_bw() +
    labs(title=ifelse(is.null(plot_title), gsub(".RData", "", data_name), plot_title),
	     x=paste0("PC1: ", round(pca_file$eig[1]/sum(pca_file$eig)*100, 2), "%"), 
         y=paste0("PC2: ", round(pca_file$eig[2]/sum(pca_file$eig)*100, 2),"%")) +
    theme(axis.title = element_text(size=7.5), axis.text = element_text(size=6), 
          axis.ticks = element_line(size=0.1),
          legend.title = element_text(size=6), legend.text = element_text(size=6), 
          legend.key.width = unit(2, "mm"),legend.key.height = unit(3, "mm"),
          legend.margin=margin(0,0,0,-8),
          title = element_blank(),
          panel.grid = element_line(size=0.05), panel.border = element_rect(size=0.2))
  
  if(print == T) {print(myplot)} else {return(myplot)}
}
 
myshape = c("pop1"=21, "pop2"=24, "pop3"=22)
datasets = list.files("./PCAdata")

pdf("all_PCAplots.pdf")
for (m in 1:length(datasets)) {
  load(paste0("./PCAdata/", datasets[m]))
  # loaded data: pca, scores, snpmiss
  plot_PCA_scores(datasets[m], scores, pca, myshape)
  rm(pca, scores, snpmiss)
}
dev.off()
rm(datasets)

## note: plots can be modified into panels using functions ggarrange() and annotate_figure(). 



######### supplementary: PCA without centering ############
glPCA_center_scale = function(SNP_data, loci_number, 
                              do_center=T, do_scale=T, ...) {
  
  gl = new("genlight", SNP_data, ploidy=2)
  pca = glPca(gl, nf=loci_number, # nloci > nindv, keep all axes
              center = do_center, scale = do_scale) 
  scores = as.data.frame(pca$scores)
  scores$indv = rownames(scores)
  scores$pop = gsub(".indv.*", "", scores$indv)
  
  indvmiss = NULL
  for (i in 1:nrow(SNP_data)) {
    x = SNP_data[i, ]
    y = length(x[x == "NA"]) / loci_number
    indvmiss = c(indvmiss, y)
  }
  scores$indvmiss = round(indvmiss, 2)
  
  snpmiss = NULL
  for (i in 1:loci_number) {
    x = SNP_data[, i]
    y = length(x[x == "NA"]) / nrow(SNP_data)
    snpmiss = c(snpmiss, y)
  }
  
  pca <<- pca
  scores <<- scores
  snpmiss <<- snpmiss
  
}

nloci = 5000 
nindv = 25

## raw SNP matrices
load("./simulation/cline_mig50_SNP.RData")
glPCA_center_scale(SNP, nloci, do_center = F)
plot = plot_PCA_scores(cline, scores, pca, myshape, plot_title="cline_raw_no-center")
save(pca, scores, snpmiss, plot, file="./PCA-nocenter/cline_mig50_PCA_nocenter.RData")
rm(pca, scores, snpmiss, plot)

load("./simulation/p3_mig50_SNP.RData")
glPCA_center_scale(SNP, nloci, do_center = F)
plot = plot_PCA_scores(p3_mig, scores, pca, myshape, plot_title="p3_mig, no centering")
save(pca, scores, snpmiss, plot, file="./PCA-nocenter/p3_mig50_PCA_nocenter.RData")
rm(pca, scores, snpmiss, plot)

## incomplete SNP matrices
SNP_data = c("cline_mig50_SNP_rand_miss0.2.RData",
             "cline_mig50_SNP_biasINDV_miss0.2.RData", 
             "cline_mig50_SNP_biasPOP_miss0.2.RData",
             "cline_mig50_SNP_biasPOP_pop3_miss0.2.RData")

for (m in 1:length(SNP_data)) {
  load(paste0("./MISSdata/", SNP_data[m]))
  
  if (m==1) {
    glPCA_center_scale(rand, nloci, do_center = F)
  }
  if (m==2) {
    glPCA_center_scale(biasINDV, nloci, do_center = F)
  }
  if (m > 2) {
    glPCA_center_scale(biasPOP, nloci, do_center = F)
  }
  
  plot = plot_PCA_scores(SNP_data[m], scores, pca, myshape)
  
  save(pca, scores, snpmiss, plot, 
       file=paste0("./PCA-nocenter/", SNP_data[m], "_PCA_nocenter.RData"))
  rm(pca, scores, snpmiss, plot)
  rm(biasINDV, biasPOP)
}

#Warning message:
#  In glDotProd(x, center = center, scale = scale, alleleAsUnit = alleleAsUnit,  :
#                 Null variances have been detected; corresponding alleles won't be standardized.





############## PCA on the empirical data sets #################
sif = read.csv("sif_ref72.csv") # samples in TableS1
sif = sif[order(sif$well),] 
# so that vcf sample IDs are in the same order as sif sample IDs
datasets = c("bat20_ref_72", "bat10_ref_72i85", "bat1_ref_72i98") 

library(vcfR) #vcfR 1.11.0 
library(adegenet) #adegenet 2.1.3 
for (i in datasets) {
  
  vcf = read.vcfR(paste0(i,".recode.vcf"))
  gl = vcfR2genlight(vcf)
  ploidy(gl) = 2
  pca = glPca(gl, nf=71) # retain all axes
  scores = as.data.frame(pca$scores) 
  
  if (all(rownames(scores) == sif$sample)) { #row names of the dataframe are sample ID, double check if in the same order as in sif
    scores$state = sif$state.abr
    scores$pop = sif$pop
    scores$indvmiss = sif[, paste0(i,".imiss")]
  }
  
  save(gl, pca, scores, file=paste0(i, "_glpca.RData"))
}

## note: empirical data are plotted using the same functions as above





