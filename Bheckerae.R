# Bathymodiolus heckerae SNP analyses

#Packages to install and access, only need to run once
install.packages('adegenet')
install.packages('hierfstat')
install.packages("factoextra")
install.packages('poppr')
install.packages('ape')
install.packages('devtools')
install.packages('dplyr')
install.packages('conStruct')
install.packages('pegas')

library('adegenet')
library('hierfstat')
library('factoextra')
library('poppr')
library('ape')
library('dplyr')
library('devtools')
library('conStruct')
library('pegas')


##Navigate to working directory where iPyrad output files are
#setwd("/path/to/output")

#Import diploid SNP data as genind object (adegenet)
#if population information is missing (NULL), it is seeked in x$pop.
#If NULL again, all individuals are assumed from the same population

##Add site information after using genind_obj@pop, code directly below 
bheckerae_data <- read.structure(file = "data2_subdata_noref.str", n.ind = 87, n.loc = 4142, onerowperind = FALSE, col.lab = 1, ask = FALSE)

##add site info, use for collection/mussel pot specific patterns
##1: mussel pot B6, 2: ROV_G, 3: SBlue_02, 4: mussel pot B1, 5: mussel pot B2, 6: mussel pot B4
sites <- as.vector(c(rep(1,33), rep(2, 1), rep(3, 1), rep(4, 17), rep(5, 10), rep(6, 25))) ##creates a vector with coded sites 1-3
population.factor <- as.factor(sites)
bheckerae_data@pop <- population.factor

#If desired to filter further, this function tests, for a series of loci, the hypothesis that genotype frequencies follow the Hardy–Weinberg equilibrium
#Uses pegas package, replaces HWE.test.genind in the adegenet package.
#hw_test <- hw.test(bchildressi_data,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix"))
#write.csv(hw_test, "HWE_results.csv")

#Convert genind object to hierfstat data frame - with population data (hierfstat)
bheckerae_pops <- genind2hierfstat(bheckerae_data)

#Run basic stats and genetic distance (hierfstat)
stats <- basic.stats(bheckerae_pops)

#overall stats:
# Ho: mean observed heterozygosities, Hs: mean gene diversities within population,
# Ht: Gene diversities overall and corrected Htp, Dst: amount of gene diversity among samples, corrected Dstp.
# Dest: a measure of population differentiation 

##If you have population data... can also run....
dist <- genet.dist(bheckerae_pops)
#Pops: 1: mussel pot B6, 2: ROV_G, 3: SBlue_02, 4: mussel pot B1, 5: mussel pot B2, 6: mussel pot B4

#Compute and save pairwise FST (hierfstat)
pairwise_results <- pairwise.neifst(bheckerae_pops, diploid=TRUE)
write.csv(pairwise_results, "pairwise_fst_results.csv")

#Principle Components Analysis (ade4)
bheckerae_data2 <- tab(bheckerae_data, NA.method="mean")
bheckerae_pca <- dudi.pca(bheckerae_data2, scale = FALSE, center= TRUE, scannf = FALSE, nf = 3)

#Visualize PCA (ade4, factoextra)
fviz_eig(bheckerae_pca)
fviz_pca_ind(bheckerae_pca, col.ind = bheckerae_data$pop, repel = TRUE, label="none")
fviz_pca_ind(bheckerae_pca, axes = c(1,2), col.ind = bheckerae_data$pop, repel = TRUE, label="none")

##############
##redo analyses on unlinked SNPs only

#Import diploid SNP data as genind object (adegenet)
#if population information is missing (NULL), it is seeked in x$pop.

##Add site information after using genind_obj@pop, code directly below 
##to determine number of fields in file,
#cat data2_subdata.ustr | awk '{ print NF}'
##number of loci = above result - 1 (because of sample ID column)
bheckerae_data <- read.structure(file = "data2_subdata_unlinked.str", n.ind = 87, n.loc = 4114, onerowperind = FALSE, col.lab = 1, ask = FALSE)

##add site info, use for collection/mussel pot specific patterns
##1: mussel pot B6, 2: ROV_G, 3: SBlue_02, 4: mussel pot B1, 5: mussel pot B2, 6: mussel pot B4
sites <- as.vector(c(rep(1,33), rep(2, 1), rep(3, 1), rep(4, 17), rep(5, 10), rep(6, 25))) ##creates a vector with coded sites 1-3
population.factor <- as.factor(sites)
bheckerae_data@pop <- population.factor

#Convert genind object to hierfstat data frame - with population data (hierfstat)
bheckerae_pops <- genind2hierfstat(bheckerae_data)

#Run basic stats and genetic distance (hierfstat)
stats <- basic.stats(bheckerae_pops)

#overall stats:
# Ho: mean observed heterozygosities, Hs: mean gene diversities within population,
# Ht: Gene diversities overall and corrected Htp, Dst: amount of gene diversity among samples, corrected Dstp.
# Dest: a measure of population differentiation 

##If you have population data... can also run....
dist <- genet.dist(bheckerae_pops)
#Pops: 1: mussel pot B6, 2: ROV_G, 3: SBlue_02, 4: mussel pot B1, 5: mussel pot B2, 6: mussel pot B4

##Compute and save pairwise FST (hierfstat)
pairwise_results <- pairwise.neifst(bheckerae_pops, diploid=TRUE)
write.csv(pairwise_results, "pairwise_fst_results.csv")

#Principle Components Analysis (ade4)
bheckerae_data2 <- tab(bheckerae_data, NA.method="mean")
bheckerae_pca <- dudi.pca(bheckerae_data2, scale = FALSE, center= TRUE, scannf = FALSE, nf = 3)

#Visualize PCA (ade4, factoextra)
fviz_eig(bheckerae_pca)
fviz_pca_ind(bheckerae_pca, col.ind = bheckerae_data$pop, repel = TRUE, label="none")
fviz_pca_ind(bheckerae_pca, axes = c(1,2), col.ind = bheckerae_data$pop, repel = TRUE, label="none")



