## G. childressi SNP analyses

## Packages to install and access, only need to run once
#install.packages('adegenet')
#install.packages('hierfstat')
#install.packages("factoextra")
#install.packages('poppr')
#install.packages('ape')
#install.packages('devtools')
#install.packages('dplyr')
#install.packages('conStruct')
#install.packages('pegas')


##load packages
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
bchildressi_data <- read.structure(file = "wgenome_20_noref.str", n.ind = 81, n.loc = 21292, onerowperind = FALSE, col.lab = 1, ask = TRUE)

##add site info
sites <- as.vector(c(rep(1,23), rep(2, 6), rep(1, 30), rep(3, 22))) ##creates a vector with coded sites 1-3
population.factor <- as.factor(sites)
bchildressi_data@pop <- population.factor

#Convert genind object to hierfstat data frame - with population data (hierfstat)
bchildressi_pops <- genind2hierfstat(bchildressi_data)

#Run basic stats and genetic distance (hierfstat)
stats <- basic.stats(bchildressi_pops)

#overall stats:
# Ho: mean observed heterozygosities, Hs: mean gene diversities within population,
# Ht: Gene diversities overall and corrected Htp, Dst: amount of gene diversity among samples, corrected Dstp.
# Dest: a measure of population differentiation 


##If you have population data... can also get genetic distance (Dch) among populations
dist <- genet.dist(bchildressi_pops)
#Pops/sites considered: 1: Norfolk canyon; 2: Chincoteague seep; 3: Baltimore canyon


#Compute and save pairwise FST (hierfstat)
pairwise_results <- pairwise.neifst(bchildressi_pops, diploid=TRUE)
write.csv(pairwise_results, "pairwise_fst_results.csv")

#Principle Components Analysis (ade4)
bchildressi_data2 <- tab(bchildressi_data, NA.method="mean")
bchildressi_pca <- dudi.pca(bchildressi_data2, scale = FALSE, center= TRUE, scannf = FALSE, nf = 3)

#Visualize PCA (ade4, factoextra)
fviz_eig(bchildressi_pca)
fviz_pca_ind(bchildressi_pca, col.ind = bchildressi_data$pop, repel = TRUE, label="ind")
fviz_pca_ind(bchildressi_pca, axes = c(1,2), col.ind = bchildressi_data$pop, repel = TRUE, label="none")


#This function tests, for a series of loci, the hypothesis that genotype frequencies follow the Hardy–Weinberg equilibrium
#Uses pegas package, replaces HWE.test.genind in the adegenet package.
hw_test <- hw.test(bchildressi_data,pop=bchildressi_data@pop,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix"))
write.csv(hw_test, "HWE_results.csv")

#check what the HWE statistic for each population is, first separate the populations with the function seppop(). 
# To focus on the analytical p-value only, set  B = 0.
(bchildressi.pop <- seppop(bchildressi_data) %>% lapply(hw.test, B = 0))

## Now we have one matrix per pop, but all we care about are the p-values, which are in the third column.
## We can use the functions sapply and [ to loop to create a matrix that only contains populations in columns and loci in rows.
(bchildressi.mat <- sapply(bchildressi.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows                                                                                                                                                

## visualizing this as a heatmap #Too many loci, doesn't work to visualize this way
## only care whether or not a given locus is in or out of HWE, we will thus define an α value and set everything above that value to 1 
## so that we can visually detect candidate loci where we might reject the Ho of HWE.
alpha  <- 0.05
newmat <- bchildressi.mat
newmat[newmat > alpha] <- 1


write.csv(newmat, "HWE_results_transformed.csv")
##used excel to seperate neutral from non-neutral loci, based on pvalues < or = 1 for each population in above transformed file (with each column representing a population)
neutral <- read.csv("neutral_loci.csv") ##read neutral loci back in AS LIST
neutrallist <- as.matrix(neutral) ##turns loci into matrix format to use for subsetting only neutral loci below


## Visualize missing data using poppr
## The poppr function info_table will help you visualize missing data so that you can assess how to treat these further 
## using the function missingno.
info_table(bchildressi_data, plot = FALSE)

## Remove loci with missing data for Bayesass analysis 
## function missingno to remove loci or individuals, 
## or replace missing data with zeroes or the average values of the locus.
## When removing loci or genotypes, you can specify a 'cutoff' representing the percent missing to be removed. The default is 0.05 (5%).
#miss <- bchildressi_data %>% missingno("loci", cutoff= 0.5) %>% info_table(plot = TRUE)
loci <- missingno(bchildressi_data, type = "loci", cutoff = 0, quiet = FALSE, freq = FALSE)
loci

##add pop/site info if not already added
sites <- as.vector(c(rep(1,23), rep(2, 6), rep(1, 30), rep(3, 22))) ##creates a vector with coded sites 1-3
population.factor <- as.factor(sites)
loci@pop <- population.factor
loci


###############
## Use R function genind2structure to write the filtered genind obj back into a structure file
source('~/Dropbox/Smithsonian/Deepsearch/bchildressi/bchildressi_toBplatifrons_outfiles_20missingonly/R_genetics_conv/genind2structure.R', encoding = 'UTF-8')
genind2structure(loci, file="wgenome_20_unlinked_filtered_structure.txt", pops=TRUE)
loci 
## The structure file then needs to be converted into immanc format for the Bayesass analysis (BA3 program) using the script pyradStr2immanc.pl
## Manually remove header row first
## wget https://github.com/stevemussmann/file_converters/blob/master/pyradStr2immanc.pl
## In command line, where popfile.txt is a tabdeliminated file of "samplename popname"
#file_converters/pyradStr2immanc.pl -m popfile.txt -s wgenome_20_unlinked_filtered_structure.txt -o wgenome_20_unlinked.immanc

## NOTE:
## checked loci number in file using the following awk command
## also check pop designation are present for every sample before moving forward
## awk '{print $3}' wgenome_20.immanc | sort | uniq | wc -l

## Then ran Bayesass analysis in command line:
## ran locally with  BA3 program (brannala), 
## need to check trace files to see if enough iterations (used Tracer program locally downloaded)
#/Users/Odontodactylus/Programs/BA3/BA3SNP -v -o wgenome_20output.txt -i 1250000 -b 150000 -g wgenome_20.immanc

####################
## in R:
#Convert filtered genind object (to remove loci with missing data) to hierfstat data frame - with population data (hierfstat) to redo stats
bchildressi_filtered <- genind2hierfstat(loci)

#Run basic stats and genetic distance (hierfstat)
stats <- basic.stats(bchildressi_filtered)
stats

## use R function above to write genind object to stru file 
## will need to convert this filtered structure file into immanc format to run Bayesass analysis, using pyradStr2immanc.pl script (or PDGSpider)
## see https://github.com/stevemussmann/file_converters for more info on conversion script
## see https://github.com/stevemussmann/BayesAss3-SNPs for more information on modified BA3 program to run the Bayesass analysis on many loci!


############ Analyses for loci under selection ###########
## SEE: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

## Convert vcf format genotype data from iPYRAD output using PLINK to be readable, ran locally 
#/Users/Danielle/Programs/plink --vcf bathy_spp_toBplatifrons.vcf --allow-extra-chr --recode --make-bed --out bathy_toBplatifrons.vcf_plink.raw

## R:
## pcadapt: statistical method implemented in pcadapt assumes that markers excessively related with population structure 
## are candidates for local adaptation; population structure determined through PCA

#install.packages("pcadapt")
library(pcadapt)

## Convert vcf format genotype data from iPYRAD output using PLINK to be readable 
#/Users/Danielle/Programs/plink --vcf bchildressi_toBplatifrons.vcf --allow-extra-chr --recode --make-bed --out bchildressi_toBplatifrons.vcf_plink.raw

##set path to PLINK coverted files
path_to_file <- "/path/to/bchildressi_toBplatifrons.vcf_plink.raw.bed"
filename <- read.pcadapt(path_to_file, type = "bed")

## First run the PCA
## To choose K, principal component analysis should first be performed with a large enough number of principal components (e.g. K=20).
x <- pcadapt(input = filename, K = 20)

#plot(x, option = "screeplot")
#plot(x, option = "screeplot", K = 7)

## steep curve cutoff shows 2 principal components, k=3

# Looking at population structure beyond K = 2 confirms the results of the scree plot. 
poplist.names <- c(rep("Norfolk",23), rep("Chincoteague", 6), rep("Norfolk", 30), rep("Baltimore", 22))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)

#The third and the fourth principal components do not ascertain population structure anymore.
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

## Computing the test statistic based on PCA
x <- pcadapt(filename, K = 3)
plot(x, option = "scores", pop = poplist.names)

#Exploring SNPs and their distribution
summary(x)
plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange") #the excess of small p-values indicates the presence of outliers.
plot(x, option = "stat.distribution") #presence of outliers is also visible when plotting a histogram of the test statistic 𝐷𝑗.

##Choosing cutoff for outlier detection

##Benjamini-Hochberg Procedure - medium conservative
#padj <- p.adjust(x$pvalues,method="BH")
#alpha <- 0.05
#outliers <- which(padj < alpha)
#length(outliers)
## 53,805 outliers, too many

## Bonferroni correction - conservative
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
## USED THIS: Conservative outlier detection = 3,429 outliers

##Linkage Disequilibrium (LD) thinning ##not an issue for RADseq data

##Association between PCs and outliers
## associate outliers with one of the K principal component to have indication about evolutionary pressure.
snp_pc <- get.pc(x, outliers)
print(snp_pc[,2])
##most seem to be asscocited with PC1?



######### Re-run PCA on unlinked SNPs only for Bayesass ############### 

#Import diploid SNP data as genind object (adegenet)
##Add site information after using genind_obj@pop

##Add site information after using genind_obj@pop, code directly below 
bchildressi_data <- read.structure(file = "wgenome_20_unlinked.str", n.ind = 81, n.loc = 21220, onerowperind = FALSE, col.lab = 1, ask = FALSE)

##add site info
sites <- as.vector(c(rep(1,23), rep(2, 6), rep(1, 30), rep(3, 22))) ##creates a vector with coded sites 1-3
population.factor <- as.factor(sites)
bchildressi_data@pop <- population.factor

#Convert genind object to hierfstat data frame - with population data (hierfstat)
bchildressi_pops <- genind2hierfstat(bchildressi_data)

#Run basic stats and genetic distance (hierfstat)
stats <- basic.stats(bchildressi_pops)

#overall stats:
# Ho: mean observed heterozygosities, Hs: mean gene diversities within population,
# Ht: Gene diversities overall and corrected Htp, Dst: amount of gene diversity among samples, corrected Dstp.
# Dest: a measure of population differentiation 

##If you have population data... can also get genetic distance (Dch) among populations
dist <- genet.dist(bchildressi_pops)
#Pops/sites considered: 1: Norfolk canyon; 2: Chincoteague seep; 3: Baltimore canyon

#Compute and save pairwise FST (hierfstat)
pairwise_results <- pairwise.neifst(bchildressi_pops, diploid=TRUE)
write.csv(pairwise_results, "pairwise_fst_results_unlinked.csv")

#Principle Components Analysis (ade4)
bchildressi_data2 <- tab(bchildressi_data, NA.method="mean")
bchildressi_pca <- dudi.pca(bchildressi_data2, scale = FALSE, center= TRUE, scannf = FALSE, nf = 3)

#Visualize PCA (ade4, factoextra)
fviz_eig(bchildressi_pca)
fviz_pca_ind(bchildressi_pca, col.ind = bchildressi_data$pop, repel = TRUE, label="ind")
fviz_pca_ind(bchildressi_pca, axes = c(1,2), col.ind = bchildressi_data$pop, repel = TRUE, label="none")


## Remove loci with too much missing data for Bayesass analysis 
## function missingno to remove loci or individuals, 
## or replace missing data with zeroes or the average values of the locus.
## When removing loci or genotypes, you can specify a 'cutoff' representing the percent missing to be removed. The default is 0.05 (5%).
loci <- missingno(bchildressi_data, type = "loci", cutoff = 0, quiet = FALSE, freq = FALSE)
loci

##add pop/site info
sites <- as.vector(c(rep(1,23), rep(2, 6), rep(1, 30), rep(3, 22))) ##creates a vector with coded sites 1-3
population.factor <- as.factor(sites)
loci@pop <- population.factor
loci

## in R:
#Convert filtered genind object (to remove loci with missing data) to hierfstat data frame - with population data (hierfstat) to redo stats
bchildressi_filtered <- genind2hierfstat(loci)

#Run basic stats and genetic distance (hierfstat)
stats <- basic.stats(bchildressi_filtered)
stats

## Then ran Bayesass analysis in command line:
## on iMAC with original BA3 program (brannala), need to check trace files to see if enough iterations
#/Users/Odontodactylus/Programs/BA3/BA3SNP -v -o wgenome_20output.txt -i 1250000 -b 250000 -g wgenome_20.immanc

#see https://rpubs.com/lbenestan/population_structure for more plot ideas

###############
##create a subset of loci based on HWE tests above
neutralloci <- loci[, loc=neutrallist] ##subsets the loci (after filtering for missing data ) to include only neutral loci
neutralloci 
##get the same results if you subset loci first, then remove missing data
#neutralloci <- bchildressi_data[, loc=neutrallist] ##subsets the loci (after filtering for missing data ) to include only neutral loci
#neutralloci
#neutralloci <- missingno(neutralloci, type = "loci", cutoff = 0.5, quiet = FALSE, freq = FALSE)

## Use R function genind2structure to write the filtered genind obj back into a structure file
#source('~/Dropbox/Smithsonian/Deepsearch/bchildressi/bchildressi_toBplatifrons_outfiles_20missingonly/R_genetics_conv/genind2structure.R', encoding = 'UTF-8')
genind2structure(neutralloci, file="wgenome_20_unlinked_filtered_structure_neutralONLY.str", pops=TRUE)

## The structure file then needs to be converted into immanc format for the Bayesass analysis (BA3 program) using the script pyradStr2immanc.pl
## wget https://github.com/stevemussmann/file_converters/blob/master/pyradStr2immanc.pl
## manually remove file header with column labels
## In command line, where popfile.txt is a tabdeliminated file of "samplename popname"
## file_converters/pyradStr2immanc.pl -m popfile.txt -s wgenome_20_unlinked_filtered_structure_neutralONLY.str -o wgenome_20_neutralloci.immanc

#################### Re-run BayesASS with neutral loci ONLY

## checked loci number in file using the following awk command 
## awk '{print $3}' wgenome_20_neutralloci.immanc | sort | uniq | wc -l

## Then ran Bayesass analysis in command line:
## on iMAC with original BA3 program (brannala), running more smoothly- need to check trace files to see if enough iterations
##/Users/Odontodactylus/Programs/BA3/BA3SNP -v -o wgenome_20output_neutral.txt -i 10000000 -b 250000 -n 100 -g wgenome_20_neutralloci.immanc
#/Users/Odontodactylus/Programs/BA3/BA3SNP -v -o wgenome_20output_neutral.txt -i 10000000 -b 250000 -n 100 -g -m 0.15 -a 0.35 -f 0.02 -t wgenome_20_neutralloci.immanc

##parameters optimized to ensure acceptance rates between 20-60% for all adjustable parameters (according to Wilson and Rannala 2003), with allele (-a) frequencies at 0.3, inbreeding coefficients (-f) at 0.02, and migration rates (-m) at 0.1 (default).
##convergence examined based on trace file, using Tracer v1.7.2