#### genotype files for the bslmmm 

#First want to check the missiging values on the bslmm, at a 5% threshold there were 4180
#123 individuals 

#This script filers and formats data to run on gemma - the progam that runs the baysesian sparse linear mixed model.
# It does so a number of ways. The filtering approach that was ultimately used in the analysis was won in which 8 individuals 
# were removed. These individuals had very high missingness values. And a shared set of loci that had a missingness of alleles at less 
# than 5% for all loci. This is all because gemma will average alleles when there is a mising value. This is tolerable for small missingness,
# but it can affect analysis at higher values. 

setwd("/media/raglandlab/ExtraDrive4/DiaClass/bslmm_pairwise/filtIndiv")

CDSD_dat <- read.table("complete_CD_SD_DiaClines_filtIndv.geno", header=F, sep = ",")
col_num <- ncol(CDSD_dat)
CDSD_genos  <- CDSD_dat[,-c(1:3)]

perc_miss <- vector()
for(i in 1:nrow(CDSD_dat)){
  perc_miss[i] <- sum(is.na(CDSD_dat[i,]))/ncol(CDSD_dat)
}

sum(perc_miss < 0.05)
#[1] 4180

#Ok, so BSLMM throws out an loci that has a missingness value about 5%, which makes sense. 
#Lets set up two runs, one where we allow BSLMM to set the threshhold on individual comparisons, 
#and another where we give the same set to each comparison and see if the results are any different. We should also run wthout filtered individuals.

#To do: 
#1. make geno files from the printed genotypes in vcf to genos. 
#2. Set up a filtered individuals and unfiltered comparsion 
#3. Set up comparisons for common set of loci and unique

#------------------------------------------#
# Make 'geno' files for bslmm
#------------------------------------------#
# This involves taking the genotype files make in vcf_to_genos.R and adding position, reference allele, and alternate allele collumns
setwd("/media/raglandlab/ExtraDrive4/Clines_3/")
library(vcfR)
#Load in the list of shared loci  
jeff_snps <- read.table("./data/position_fix.txt")
jeff_snps <- as.vector(jeff_snps$V1)

#Load in the list of polarized loci 
load("./data/polarized_loci.Rdata")

#Load in the vcf to get the position, reference, and alternate that the 
Das <- read.vcfR("./data/snps.GATK.all.filtered.vcf", nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = T)
PLs <- extract.gt(Das, element = "PL", mask = F, as.numeric = F, return.alleles = F, IDtoRowNames = TRUE, extract = T)
Refs <- getREF(Das)
Alts <- getALT(Das)
Pos <- rownames(PLs)
pos_alleles <- data.frame(pos = Pos, refs = Refs, alts = Alts)
pos_alleles$pos <- as.character(pos_alleles$pos)
pos_alleles$refs <- as.character(pos_alleles$refs)
pos_alleles$alts <- as.character(pos_alleles$alts)

pos_alleles <- pos_alleles[pos_alleles$pos %in% jeff_snps, ]
#Ok, need to flip the alleles based on the polarization scheme. 
alleles_scored <- matrix(nrow = nrow(pos_alleles), ncol = 2)
for(i in 1:length(pol_ind)){
  if(pol_ind[i]==FALSE){
    alleles_scored[i,1] <- pos_alleles$refs[i]
    alleles_scored[i,2] <- pos_alleles$alts[i]
  }else if(pol_ind[i]==TRUE){
    alleles_scored[i,1] <- pos_alleles$alts[i]
    alleles_scored[i,2] <- pos_alleles$refs[i]
  }
  
}
rownames(alleles_scored) <- pos_alleles$pos
colnames(alleles_scored) <- c("Ref", "Alt")
alleles_scored <- as.data.frame(alleles_scored, stringsAsFactors = F)
#Load up the genotypes and then merge them. 
SD_dat <- read.table("./data/SD_genosAll_polarized.txt", header=T)
SD_genos <- merge(alleles_scored, SD_dat, by = "row.names")

#CD
CD_dat <- read.table("./data/CD_genosAll_polarized.txt")
CD_genos <- merge(alleles_scored, CD_dat, by = "row.names")

#ND 
ND_dat <- read.table("./data/ND_genosAll_polarized.txt")
ND_genos <- merge(alleles_scored, ND_dat, by = "row.names")

#Ok, now we'll see how many loci are maintained at missingness above below 5% in each pairwise comparison

SDCD_ind <- vector()
for(i in 1:nrow(SD_dat)){
  SDCD_ind[i] <- sum(is.na(cbind(SD_dat[i,],CD_dat[i,])))/ncol(cbind(SD_dat[i,],CD_dat[i,]))
}

SDND_ind <- vector()
for(i in 1:nrow(SD_dat)){
  SDND_ind[i] <- sum(is.na(cbind(SD_dat[i,],ND_dat[i,])))/ncol(cbind(SD_dat[i,],ND_dat[i,]))
}

CDND_ind <- vector()
for(i in 1:nrow(SD_dat)){
  CDND_ind[i] <- sum(is.na(cbind(CD_dat[i,],ND_dat[i,])))/ncol(cbind(CD_dat[i,],ND_dat[i,]))
}

sum(SDCD_ind < 0.05)
#[1] 1295, that is much lower than the ~4k that should be there. 

#Check to see if removing 8 individuals makes the difference. 
SDCD_test_mat <- merge(SD_dat, CD_dat, by="row.names")

#filter the individuals: 
rm_ind <- c("index_10nt_41", "index_10nt_47", "index_8nt_19", "index_8nt_28", "index_8nt_30", "index_8nt_32", "index_8nt_34", "index_9nt_21")
rm_x <- colnames(SDCD_test_mat) %in% rm_ind

perc_miss <- vector()
for(i in 1:nrow(SDCD_test_mat)){
  perc_miss[i] <- sum(is.na(SDCD_test_mat[i,!rm_x]))/ncol(SDCD_test_mat[,!rm_x])
}

sum(perc_miss < 0.05)
#[1] 4180, ok, so it comes down to those removed individuals. Which were removed because they had very high amounts of missingness. 

#Ok, now that that's sorted, we"ll make: 
#------------------------------------------#
# Filtered invidiviuals full set
#------------------------------------------#
# This is setting up .geno and .pheno files for bslmm using individuals filtered and full loci set, meaning bslmm will filter based on 
# the 'missingness' setting, usually 0.05 
#A genotype file for bslmm is comma separated with position, maj, min, followed by mean genotype probabilities. Each column represents a sample. 
# The genotype files contain all data to be analyzed, so all populations. 

#First we'll filter all the individuals from each group
#SD 
rm_SD <- colnames(SD_genos) %in% rm_ind
SD_genos <- SD_genos[,!rm_SD]
ncol(SD_genos) - 3
# number of individuals still present
# 64
#CD 
rm_CD <- colnames(CD_genos) %in% rm_ind
CD_genos <- CD_genos[,!rm_CD]
ncol(CD_genos) - 3
# number of individuals still present
#59
#ND 
rm_ND <- colnames(ND_genos) %in% rm_ind
ND_genos <- ND_genos[,!rm_ND]
ncol(ND_genos) - 3
# number of individuals still present 
#61

# SD vs. CD 
#genotype file 
SDCD_filtIndiv_fullSet <- merge(SD_genos, CD_genos, by = c("Row.names", "Ref", "Alt"))
#phenotype file - coded as 1 for SD and 0 for ND
SD <- rep(1, 64)
CD <- rep(0, 59)
SDCD_filtIndiv_fullSet_pheno <- matrix(c(SD, CD), ncol = 1)

write.table(SDCD_filtIndiv_fullSet, file = "./data/SDCD_filtIndiv_fullSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDCD_filtIndiv_fullSet_pheno, file = "./data/SDCD_filtIndiv_fullSet.pheno", row.names = F, col.names = F, quote = F)

# SD vs. ND 

SDND_filtIndiv_fullSet <- merge(SD_genos, ND_genos, by = c("Row.names", "Ref", "Alt"))
#phenotype file - coded as 1 for SD and 0 for ND
SD <- rep(1, 64)
ND <- rep(0, 61)
SDND_filtIndiv_fullSet_pheno <- matrix(c(SD, ND), ncol = 1)

write.table(SDND_filtIndiv_fullSet, file = "./data/SDND_filtIndiv_fullSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDND_filtIndiv_fullSet_pheno, file = "./data/SDND_filtIndiv_fullSet.pheno", row.names = F, col.names = F, quote = F)

# CD vs. ND 

CDND_filtIndiv_fullSet <- merge(CD_genos, ND_genos, by = c("Row.names", "Ref", "Alt"))
#phenotype file - coded as 1 for CD and 0 for ND
CD <- rep(1, 59)
ND <- rep(0, 61)
CDND_filtIndiv_fullSet_pheno <- matrix(c(CD, ND), ncol = 1)

write.table(CDND_filtIndiv_fullSet, file = "./data/CDND_filtIndiv_fullSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(CDND_filtIndiv_fullSet_pheno, file = "./data/CDND_filtIndiv_fullSet.pheno", row.names = F, col.names = F, quote = F)

#------------------------------------------#
# Filtered invidiviuals combined set
#------------------------------------------#

#We already have the fullset genotype files, so we just need to make indexes of them and then a final index

SDCD_ind <- vector()
for(i in 1:nrow(SD_dat)){
  SDCD_ind[i] <- sum(is.na(SDCD_filtIndiv_fullSet[i,-c(1:3)]))/ncol(SDCD_filtIndiv_fullSet[i,-c(1:3)])
}
sum(SDCD_ind < 0.05)
#[1] 4180

SDND_ind <- vector()
for(i in 1:nrow(SD_dat)){
  SDND_ind[i] <- sum(is.na(SDND_filtIndiv_fullSet[i,-c(1:3)]))/ncol(SDND_filtIndiv_fullSet[i,-c(1:3)])
}
sum(SDND_ind < 0.05)
#[1] 4123
CDND_ind <- vector()
for(i in 1:nrow(SD_dat)){
  CDND_ind[i] <- sum(is.na(CDND_filtIndiv_fullSet[i,-c(1:3)]))/ncol(CDND_filtIndiv_fullSet[i,-c(1:3)])
}
sum(CDND_ind < 0.05)
#[1] 4432

#Ok, how many are shared across the whole set
full_ind <- SDCD_ind < 0.05 & SDND_ind < 0.05 & CDND_ind < 0.05
sum(full_ind)
#[1] 3611

#don't really need to print a new phenotype file, but to keep things straight and un confusing 
write.table(SDCD_filtIndiv_fullSet[full_ind,], file = "./data/SDCD_filtIndiv_filtSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDCD_filtIndiv_fullSet_pheno, file = "./data/SDCD_filtIndiv_filtSet.pheno", row.names = F, col.names = F, quote = F)

write.table(SDND_filtIndiv_fullSet[full_ind,], file = "./data/SDND_filtIndiv_filtSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDND_filtIndiv_fullSet_pheno, file = "./data/SDND_filtIndiv_filtSet.pheno", row.names = F, col.names = F, quote = F)

write.table(CDND_filtIndiv_fullSet[full_ind,], file = "./data/CDND_filtIndiv_filtSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(CDND_filtIndiv_fullSet_pheno, file = "./data/CDND_filtIndiv_filtSet.pheno", row.names = F, col.names = F, quote = F)

#------------------------------------------#
# UnFiltered invidiviuals full set
#------------------------------------------#
#Load up the genotypes and then merge them. 
SD_dat <- read.table("./data/SD_genosAll_polarized.txt", header=T)
SD_genos <- merge(alleles_scored, SD_dat, by = "row.names")

#CD
CD_dat <- read.table("./data/CD_genosAll_polarized.txt")
CD_genos <- merge(alleles_scored, CD_dat, by = "row.names")

#ND 
ND_dat <- read.table("./data/ND_genosAll_polarized.txt")
ND_genos <- merge(alleles_scored, ND_dat, by = "row.names")

# SD vs. CD 
#genotype file 
SDCD_fullIndiv_fullSet <- merge(SD_genos, CD_genos, by = c("Row.names", "Ref", "Alt"))
#phenotype file - coded as 1 for SD and 0 for ND
SD <- rep(1, 64)
CD <- rep(0, 64)
SDCD_fullIndiv_fullSet_pheno <- matrix(c(SD, CD), ncol = 1)

write.table(SDCD_fullIndiv_fullSet, file = "./data/SDCD_fullIndiv_fullSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDCD_fullIndiv_fullSet_pheno, file = "./data/SDCD_fullIndiv_fullSet.pheno", row.names = F, col.names = F, quote = F)

# SD vs. ND 

SDND_fullIndiv_fullSet <- merge(SD_genos, ND_genos, by = c("Row.names", "Ref", "Alt"))
#phenotype file - coded as 1 for SD and 0 for ND
SD <- rep(1, 64)
ND <- rep(0, 64)
SDND_fullIndiv_fullSet_pheno <- matrix(c(SD, ND), ncol = 1)

write.table(SDND_fullIndiv_fullSet, file = "./data/SDND_fullIndiv_fullSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDND_fullIndiv_fullSet_pheno, file = "./data/SDND_fullIndiv_fullSet.pheno", row.names = F, col.names = F, quote = F)

# CD vs. ND 

CDND_fullIndiv_fullSet <- merge(CD_genos, ND_genos, by = c("Row.names", "Ref", "Alt"))
#phenotype file - coded as 1 for CD and 0 for ND
CD <- rep(1, 64)
ND <- rep(0, 64)
CDND_fullIndiv_fullSet_pheno <- matrix(c(CD, ND), ncol = 1)

write.table(CDND_fullIndiv_fullSet, file = "./data/CDND_fullIndiv_fullSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(CDND_fullIndiv_fullSet_pheno, file = "./data/CDND_fullIndiv_fullSet.pheno", row.names = F, col.names = F, quote = F)

#------------------------------------------#
# UnFiltered invidiviuals combined set
#------------------------------------------#


SDCD_ind <- vector()
for(i in 1:nrow(SD_dat)){
  SDCD_ind[i] <- sum(is.na(SDCD_fullIndiv_fullSet[i,-c(1:3)]))/ncol(SDCD_fullIndiv_fullSet[i,-c(1:3)])
}
sum(SDCD_ind < 0.05)
#[1] 4180

SDND_ind <- vector()
for(i in 1:nrow(SD_dat)){
  SDND_ind[i] <- sum(is.na(SDND_fullIndiv_fullSet[i,-c(1:3)]))/ncol(SDND_fullIndiv_fullSet[i,-c(1:3)])
}
sum(SDND_ind < 0.05)
#[1] 4123
CDND_ind <- vector()
for(i in 1:nrow(SD_dat)){
  CDND_ind[i] <- sum(is.na(CDND_fullIndiv_fullSet[i,-c(1:3)]))/ncol(CDND_fullIndiv_fullSet[i,-c(1:3)])
}
sum(CDND_ind < 0.1)
#[1] 4432

#Ok, how many are shared across the whole set, we'll raise the CDND comparison to 10% because of the handful of individuals 
# across those groups with low coverage. 
full_ind <- SDCD_ind < 0.05 & SDND_ind < 0.05 & CDND_ind < 0.1
sum(full_ind)
#[1] 1261

#don't really need to print a new phenotype file, but to keep things straight and un confusing 
write.table(SDCD_fullIndiv_fullSet[full_ind,], file = "./data/SDCD_fullIndiv_filtSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDCD_fullIndiv_fullSet_pheno, file = "./data/SDCD_fullIndiv_filtSet.pheno", row.names = F, col.names = F, quote = F)

write.table(SDND_fullIndiv_fullSet[full_ind,], file = "./data/SDND_fullIndiv_filtSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(SDND_fullIndiv_fullSet_pheno, file = "./data/SDND_fullIndiv_filtSet.pheno", row.names = F, col.names = F, quote = F)

write.table(CDND_fullIndiv_fullSet[full_ind,], file = "./data/CDND_fullIndiv_filtSet.geno", row.names = F, col.names = F, quote = F, sep = ",")
write.table(CDND_fullIndiv_fullSet_pheno, file = "./data/CDND_fullIndiv_filtSet.pheno", row.names = F, col.names = F, quote = F)
