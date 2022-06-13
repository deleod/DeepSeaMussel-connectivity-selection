## Bathymodioline species analyses for loci under selection
## SEE: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

## Convert vcf format genotype data from iPYRAD output using PLINK to be readable, run locally 
/Users/Danielle/Programs/plink --vcf bathy_spp_toBplatifrons.vcf --allow-extra-chr --recode --make-bed --out bathy_toBplatifrons.vcf_plink.raw

## R:
## pcadapt: statistical method implemented in pcadapt assumes that markers excessively related with population structure 
## are candidates for local adaptation; population structure determined through PCA

#install.packages("pcadapt") ##only need to run once
library(pcadapt)

setwd("~/path/to/converted/output/files")

#convert your genotype file to the bed format, which is PLINK binary biallelic genotype table
path_to_file <- "/Users/Danielle/Dropbox/Smithsonian/Deepsearch/bathy_spp/bathy_toBplatifrons.vcf_plink.raw.bed"
filename <- read.pcadapt(path_to_file, type = "bed")

## First run the PCA
## To choose K, principal component analysis should first be performed with a large enough number of principal components (e.g. K=20).
x <- pcadapt(input = filename, K = 20)

plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 5)

  ## steep curve cutoff shows 2 principal components, k=2

# Looking at population structure beyond K = 2 confirms the results of the scree plot. 
#poplist.names <- c(rep("Bhec", 10),rep("Bchi", 9), rep("Bhec", 1), rep("Bchi", 19), rep("Bhec", 17))
poplist.names <- c(rep("Bhec", 10),rep("Bchi_NC", 5), rep("Bchi_CT", 4), rep("Bhec", 1), rep("Bchi_NC", 9), rep("Bchi_BC", 10),rep("Bhec", 17))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)

#The third and the fourth principal components do not ascertain population structure anymore.
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

## Computing the test statistic based on PCA
x <- pcadapt(filename, K = 2)

#Exploring SNPs and their distribution
summary(x)
write.csv(x, "pcadapt_output")
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
  ## 53,805 outliers

##Bonferroni correction - conservative
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
  ## USED THIS: Conservative outlier detection = 3,429 outliers
write.csv(outliers, "pcadapt_outiers.csv")

##Association between PCs and outliers
## associate outliers with one of the K principal component to have indication about evolutionary pressure.
snp_pc <- get.pc(x, outliers)
print(snp_pc[,2])
    ##most seem to be asscoaited with PC1?
