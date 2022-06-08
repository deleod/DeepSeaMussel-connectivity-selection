## G.childressi sibling analyses

#install.packages("sequoia") ##only required once
library(sequoia) ##load package

##Navigate to appropriate working directory where ipyrad output files are (.vcf)
#setwd("/path/to/workingdirectory")

###### In terminal/locally
## Convert vcf format genotype data from iPYRAD output using PLINK to be readable for SEQUOIA (done in terminal)
/Users/Danielle/Programs/plink --vcf wgenome_20.vcf --allow-extra-chr --recode --out bchildressi_toBplatifrons.vcf_plink.raw

# Based on sequoia userguide, created subset of SNPs for analyses (done in terminal)
/Users/Danielle/Programs/plink --file bchildressi_toBplatifrons.vcf_plink.raw --geno 0.2 --maf 0.3 --indep 50 5 1 --allow-extra-chr 

## (R) SNPs with a missingness below 0.2, a minor allele frequency of at least 0.3, and which in a window of 50 SNPs, sliding by 5 SNPs per step, have a VIF of maximum 2
## advised to ‘tweak’ the parameter values until a set with a few hundred SNPs (300–700) is created

# Extract subset above (done in terminal)
/Users/Danielle/Programs/plink --file bchildressi_toBplatifrons.vcf_plink.raw --extract plink.prune.in --recodeA --allow-extra-chr --out inputfile_for_sequoia
########################

##In R:
Geno <- GenoConvert(InFile ="inputfile_for_sequoia.raw", InFormat="raw") ##subset of approx. 669 SNPs
#CheckGeno(Geno)
#
# read in lifehistory data: ID-Sex-birthyear, column names ignored
# optional: minimum & maximum birth year, when not exactly known

LH <- read.table("lifehistory.txt", header=T)
#
# duplicate check & parentage assignment (takes few minutes)

ParOUT <- sequoia(GenoM = Geno,  LifeHistData = LH,
                  MaxSibIter = 20, Err=0.005, FindMaybeRel = TRUE,
                  quiet = FALSE)
# inspect duplicates (intentional or accidental) & remove if needed
ParOUT$DupGenotype


stats <- SnpStats(Geno)

## run full pedigree reconstruction (may take up to a few hours)
SeqOUT <- sequoia(GenoM = Geno,
                  LifeHistData = LH,
                  Err = 0.001,
                  MaxSibIter = 40,
                  FindMaybeRel = TRUE, quiet = FALSE)
#
# inspect assigned parents, proportion dummy parents, etc.
SummarySeq(SeqOUT)

Relatives <- SeqOUT$MaybeRel

##determines pairwise realtionships

# save results
save(Relatives, file="Sequoia_relatives.txt")
save(SeqOUT, file="Sequoia_output_032621.RData")
writeSeq(SeqList = SeqOUT, GenoM = Geno, folder = "Sequoia-OUT")
