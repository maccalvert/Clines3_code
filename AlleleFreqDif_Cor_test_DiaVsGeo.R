# Alelle frequency difference correlations between diapause intensity and the other diapause phenotypes previously genotyped in other studies. 

#Analyses included Diapause intensity vs. Host race

#Analyses performed for each chromosome and LD group. Tables and figures printed at the end of each comparison. 

#This correlation analysis was performed using a permutation analysis. Here, permuted allele frequency difference for each "experiment" were 
#created. Correlation estimates are then calculated between permuted allele frequency differences for the two experiments under comparison. This 
#generates a null distribution of correlation coefficients to compare the empirical correlation coefficient against. 
# Because allele frequency distributions don't follow a normal distribution, we used spearman rank correlation coefficients. 

#Set up the workspace 
setwd("/media/raglandlab/ExtraDrive4/Clines_3/")

#Load in R functions
source("./src/clines3_functions.R")

library(iterators)
library(doParallel)
library(foreach)
library(ggplot2)
library(viridis)

#Load in the position information index 
load("./data/Mapped_RAD_loci.Rdata")

################################################################
# 3. Diapause intensity vs. Geography (Grant, MI - Urbana, IL) #
################################################################

#Hawthorn host race 

#Load in the empirical allele frequency differences for each 
DiaInt <- read.table("./results/SDvCDND_freqDifs.txt", header=T)
HawGeo <- read.table("./results/Haw_geo_freqDifs.txt", header=T)

#Make empty matrices to store the p_values and corelation estimates 
DiavHawGeo_pvals <- matrix(nrow=4, ncol =6)
rownames(DiavHawGeo_pvals) <- c("All_Snps", "High_LD", "Int. LD", "Low LD")
colnames(DiavHawGeo_pvals) <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_all")
DiavHawGeo_cors <- matrix(nrow=4, ncol =6)
rownames(DiavHawGeo_cors) <- c("All_Snps", "High_LD", "Int. LD", "Low LD")
colnames(DiavHawGeo_cors) <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_all")

#Chromosome 1 

SD_mat <- read.table("./data/SD_genosChr1_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr1_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr1_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_ind_p], HawGeo$freqDif[chr1_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[1,1] <- pval

DiavHawGeo_cors[1,1] <- corEst$estimate

#Chromosome 1 High LD

SD_mat <- read.table("./data/SD_genosChr1_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDH_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr1_LDH_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr1_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_H_ind_p], HawGeo$freqDif[chr1_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[2,1] <- pval

DiavHawGeo_cors[2,1] <- corEst$estimate

#Chromosome 1 int LD

SD_mat <- read.table("./data/SD_genosChr1_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDM_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr1_LDM_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr1_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_M_ind_p], HawGeo$freqDif[chr1_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[3,1] <- pval

DiavHawGeo_cors[3,1] <- corEst$estimate


#Chromosome 1 low LD

SD_mat <- read.table("./data/SD_genosChr1_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDL_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr1_LDL_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr1_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_L_ind_p], HawGeo$freqDif[chr1_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[4,1] <- pval

DiavHawGeo_cors[4,1] <- corEst$estimate


#Chromosome 2

SD_mat <- read.table("./data/SD_genosChr2_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr2_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr2_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_ind_p], HawGeo$freqDif[chr2_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[1,2] <- pval

DiavHawGeo_cors[1,2] <- corEst$estimate

#Chromosome 2 High LD

SD_mat <- read.table("./data/SD_genosChr2_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDH_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr2_LDH_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr2_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_H_ind_p], HawGeo$freqDif[chr2_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[2,2] <- pval

DiavHawGeo_cors[2,2] <- corEst$estimate

#Chromosome 2 int LD

SD_mat <- read.table("./data/SD_genosChr2_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDM_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr2_LDM_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr2_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_M_ind_p], HawGeo$freqDif[chr2_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[3,2] <- pval

DiavHawGeo_cors[3,2] <- corEst$estimate


#Chromosome 2 low LD

SD_mat <- read.table("./data/SD_genosChr2_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDL_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr2_LDL_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr2_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_L_ind_p], HawGeo$freqDif[chr2_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[4,2] <- pval

DiavHawGeo_cors[4,2] <- corEst$estimate

#Chromosome 3

SD_mat <- read.table("./data/SD_genosChr3_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr3_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr3_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_ind_p], HawGeo$freqDif[chr3_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[1,3] <- pval

DiavHawGeo_cors[1,3] <- corEst$estimate

#Chromosome 3 High LD

SD_mat <- read.table("./data/SD_genosChr3_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDH_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr3_LDH_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr3_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_H_ind_p], HawGeo$freqDif[chr3_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[2,3] <- pval

DiavHawGeo_cors[2,3] <- corEst$estimate

#Chromosome 3 int LD

SD_mat <- read.table("./data/SD_genosChr3_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDM_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr3_LDM_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr3_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_M_ind_p], HawGeo$freqDif[chr3_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[3,3] <- pval

DiavHawGeo_cors[3,3] <- corEst$estimate


#Chromosome 3 low LD

SD_mat <- read.table("./data/SD_genosChr3_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDL_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr3_LDL_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr3_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_L_ind_p], HawGeo$freqDif[chr3_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[4,3] <- pval

DiavHawGeo_cors[4,3] <- corEst$estimate

#Chromosome 4

SD_mat <- read.table("./data/SD_genosChr4_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr4_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr4_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_ind_p], HawGeo$freqDif[chr4_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[1,4] <- pval

DiavHawGeo_cors[1,4] <- corEst$estimate

#Chromosome 4 High LD

SD_mat <- read.table("./data/SD_genosChr4_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDH_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr4_LDH_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr4_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_H_ind_p], HawGeo$freqDif[chr4_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[2,4] <- pval

DiavHawGeo_cors[2,4] <- corEst$estimate

#Chromosome 4 int LD

SD_mat <- read.table("./data/SD_genosChr4_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDM_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr4_LDM_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr4_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_M_ind_p], HawGeo$freqDif[chr4_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[3,4] <- pval

DiavHawGeo_cors[3,4] <- corEst$estimate


#Chromosome 4 low LD

SD_mat <- read.table("./data/SD_genosChr4_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDL_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr4_LDL_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr4_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_L_ind_p], HawGeo$freqDif[chr4_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[4,4] <- pval

DiavHawGeo_cors[4,4] <- corEst$estimate

#Chromosome 5

SD_mat <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr5_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr5_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_ind_p], HawGeo$freqDif[chr5_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[1,5] <- pval

DiavHawGeo_cors[1,5] <- corEst$estimate

#Chromosome 5 High LD

SD_mat <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr5_LDH_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr5_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_H_ind_p], HawGeo$freqDif[chr5_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[2,5] <- pval

DiavHawGeo_cors[2,5] <- corEst$estimate

#Chromosome 5 int LD

SD_mat <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr5_LDM_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr5_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_M_ind_p], HawGeo$freqDif[chr5_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[3,5] <- pval

DiavHawGeo_cors[3,5] <- corEst$estimate


#Chromosome 5 low LD

SD_mat <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genosChr5_LDL_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genosChr5_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_L_ind_p], HawGeo$freqDif[chr5_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[4,5] <- pval

DiavHawGeo_cors[4,5] <- corEst$estimate


#All Mapped loci 
SD_mat <- read.table("./data/SD_genosAllChr_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosAllChr_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genos_AllChr_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genos_AllChr_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_chr_ind_p], HawGeo$freqDif[all_chr_ind_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[1,6] <- pval

DiavHawGeo_cors[1,6] <- corEst$estimate

#All high LD 
SD_mat <- read.table("./data/SD_genos_AllLDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDH_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genos_AllLDH_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genos_AllLDH_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDH_p], HawGeo$freqDif[all_LDH_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[2,6] <- pval

DiavHawGeo_cors[2,6] <- corEst$estimate

#All int. LD 
SD_mat <- read.table("./data/SD_genos_AllLDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDM_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genos_AllLDM_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genos_AllLDM_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDM_p], HawGeo$freqDif[all_LDM_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[3,6] <- pval

DiavHawGeo_cors[3,6] <- corEst$estimate


#All low LD
SD_mat <- read.table("./data/SD_genos_AllLDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDL_polarized.txt", header=T)
haw_7 <- read.table("./data/haw7_genos_AllLDL_polarized.txt", header=T)
UrbHaw <- read.table("./data/UrbHaw_genos_AllLDL_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDL_p], HawGeo$freqDif[all_LDL_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = haw_7, ecl2 = UrbHaw, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHawGeo_pvals[4,6] <- pval

DiavHawGeo_cors[4,6] <- corEst$estimate

# Write the estimate and pvalue files
write.table(DiavHawGeo_cors, file = "./results/DiavHawGeo_corEst.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(DiavHawGeo_pvals, file = "./results/DiavHawGeo_perPvals.txt", row.names = T, col.names = T, quote = F, sep = "\t")
#######################################
# Figures
#######################################

#Figure 5 and supplementary figure 3

#Chromosome 1 

dat <- as.data.frame(cbind(HawGeo$freqDif, DiaInt$freqDif))
lowld <- as.data.frame(cbind(HawGeo$freqDif[chr1_L_ind_p], DiaInt$freqDif[chr1_L_ind_p]))
midld <- as.data.frame(cbind(HawGeo$freqDif[chr1_M_ind_p], DiaInt$freqDif[chr1_M_ind_p]))
highld <- as.data.frame(cbind(HawGeo$freqDif[chr1_H_ind_p], DiaInt$freqDif[chr1_H_ind_p]))
colnames(dat) <- c("HawGeo", "Dia")
colnames(lowld) <- c("HawGeo", "Dia")
colnames(midld) <- c("HawGeo", "Dia")
colnames(highld) <- c("HawGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHawGeo_chr1.tiff")
p <- ggplot(dat, aes(x=HawGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HawGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HawGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HawGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 2 
lowld <- as.data.frame(cbind(HawGeo$freqDif[chr2_L_ind_p], DiaInt$freqDif[chr2_L_ind_p]))
midld <- as.data.frame(cbind(HawGeo$freqDif[chr2_M_ind_p], DiaInt$freqDif[chr2_M_ind_p]))
highld <- as.data.frame(cbind(HawGeo$freqDif[chr2_H_ind_p], DiaInt$freqDif[chr2_H_ind_p]))
colnames(dat) <- c("HawGeo", "Dia")
colnames(lowld) <- c("HawGeo", "Dia")
colnames(midld) <- c("HawGeo", "Dia")
colnames(highld) <- c("HawGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHawGeo_chr2.tiff")
p <- ggplot(dat, aes(x=HawGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HawGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HawGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HawGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 3 
lowld <- as.data.frame(cbind(HawGeo$freqDif[chr3_L_ind_p], DiaInt$freqDif[chr3_L_ind_p]))
midld <- as.data.frame(cbind(HawGeo$freqDif[chr3_M_ind_p], DiaInt$freqDif[chr3_M_ind_p]))
highld <- as.data.frame(cbind(HawGeo$freqDif[chr3_H_ind_p], DiaInt$freqDif[chr3_H_ind_p]))
colnames(dat) <- c("HawGeo", "Dia")
colnames(lowld) <- c("HawGeo", "Dia")
colnames(midld) <- c("HawGeo", "Dia")
colnames(highld) <- c("HawGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHawGeo_chr3.tiff")
p <- ggplot(dat, aes(x=HawGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HawGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HawGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HawGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 4
lowld <- as.data.frame(cbind(HawGeo$freqDif[chr4_L_ind_p], DiaInt$freqDif[chr4_L_ind_p]))
midld <- as.data.frame(cbind(HawGeo$freqDif[chr4_M_ind_p], DiaInt$freqDif[chr4_M_ind_p]))
highld <- as.data.frame(cbind(HawGeo$freqDif[chr4_H_ind_p], DiaInt$freqDif[chr4_H_ind_p]))
colnames(dat) <- c("HawGeo", "Dia")
colnames(lowld) <- c("HawGeo", "Dia")
colnames(midld) <- c("HawGeo", "Dia")
colnames(highld) <- c("HawGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHawGeo_chr4.tiff")
p <- ggplot(dat, aes(x=HawGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HawGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HawGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HawGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 5
lowld <- as.data.frame(cbind(HawGeo$freqDif[chr5_L_ind_p], DiaInt$freqDif[chr5_L_ind_p]))
midld <- as.data.frame(cbind(HawGeo$freqDif[chr5_M_ind_p], DiaInt$freqDif[chr5_M_ind_p]))
highld <- as.data.frame(cbind(HawGeo$freqDif[chr5_H_ind_p], DiaInt$freqDif[chr5_H_ind_p]))
colnames(dat) <- c("HawGeo", "Dia")
colnames(lowld) <- c("HawGeo", "Dia")
colnames(midld) <- c("HawGeo", "Dia")
colnames(highld) <- c("HawGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHawGeo_chr5.tiff")
p <- ggplot(dat, aes(x=HawGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HawGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HawGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HawGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()


#####################################################
# Diapause intensity vs. Apple Pre winter selection #
#####################################################

#Load in the empirical allele frequency differences for each 
DiaInt <- read.table("./results/SDvCDND_freqDifs.txt", header=T)
AppleGeo <- read.table("./results/Apple_geo_freqDifs.txt", header=T)

#Make empty matrices to store the p_values and corelation estimates 
DiavAppleGeo_pvals <- matrix(nrow=4, ncol =6)
rownames(DiavAppleGeo_pvals) <- c("All_Snps", "High_LD", "Int. LD", "Low LD")
colnames(DiavAppleGeo_pvals) <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_all")
DiavAppleGeo_cors <- matrix(nrow=4, ncol =6)
rownames(DiavAppleGeo_cors) <- c("All_Snps", "High_LD", "Int. LD", "Low LD")
colnames(DiavAppleGeo_cors) <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_all")

#Chromosome 1 

SD_mat <- read.table("./data/SD_genosChr1_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr1_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr1_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_ind_p], AppleGeo$freqDif[chr1_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[1,1] <- pval

DiavAppleGeo_cors[1,1] <- corEst$estimate

#Chromosome 1 High LD

SD_mat <- read.table("./data/SD_genosChr1_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDH_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr1_LDH_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr1_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_H_ind_p], AppleGeo$freqDif[chr1_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[2,1] <- pval

DiavAppleGeo_cors[2,1] <- corEst$estimate

#Chromosome 1 int LD

SD_mat <- read.table("./data/SD_genosChr1_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDM_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr1_LDM_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr1_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_M_ind_p], AppleGeo$freqDif[chr1_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[3,1] <- pval

DiavAppleGeo_cors[3,1] <- corEst$estimate


#Chromosome 1 low LD

SD_mat <- read.table("./data/SD_genosChr1_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDL_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr1_LDL_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr1_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_L_ind_p], AppleGeo$freqDif[chr1_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[4,1] <- pval

DiavAppleGeo_cors[4,1] <- corEst$estimate


#Chromosome 2

SD_mat <- read.table("./data/SD_genosChr2_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr2_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr2_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_ind_p], AppleGeo$freqDif[chr2_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[1,2] <- pval

DiavAppleGeo_cors[1,2] <- corEst$estimate

#Chromosome 2 High LD

SD_mat <- read.table("./data/SD_genosChr2_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDH_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr2_LDH_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr2_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_H_ind_p], AppleGeo$freqDif[chr2_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[2,2] <- pval

DiavAppleGeo_cors[2,2] <- corEst$estimate

#Chromosome 2 int LD

SD_mat <- read.table("./data/SD_genosChr2_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDM_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr2_LDM_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr2_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_M_ind_p], AppleGeo$freqDif[chr2_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[3,2] <- pval

DiavAppleGeo_cors[3,2] <- corEst$estimate


#Chromosome 2 low LD

SD_mat <- read.table("./data/SD_genosChr2_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDL_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr2_LDL_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr2_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_L_ind_p], AppleGeo$freqDif[chr2_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[4,2] <- pval

DiavAppleGeo_cors[4,2] <- corEst$estimate

#Chromosome 3

SD_mat <- read.table("./data/SD_genosChr3_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr3_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr3_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_ind_p], AppleGeo$freqDif[chr3_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[1,3] <- pval

DiavAppleGeo_cors[1,3] <- corEst$estimate

#Chromosome 3 High LD

SD_mat <- read.table("./data/SD_genosChr3_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDH_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr3_LDH_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr3_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_H_ind_p], AppleGeo$freqDif[chr3_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[2,3] <- pval

DiavAppleGeo_cors[2,3] <- corEst$estimate

#Chromosome 3 int LD

SD_mat <- read.table("./data/SD_genosChr3_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDM_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr3_LDM_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr3_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_M_ind_p], AppleGeo$freqDif[chr3_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[3,3] <- pval

DiavAppleGeo_cors[3,3] <- corEst$estimate


#Chromosome 3 low LD

SD_mat <- read.table("./data/SD_genosChr3_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDL_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr3_LDL_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr3_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_L_ind_p], AppleGeo$freqDif[chr3_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[4,3] <- pval

DiavAppleGeo_cors[4,3] <- corEst$estimate

#Chromosome 4

SD_mat <- read.table("./data/SD_genosChr4_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr4_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr4_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_ind_p], AppleGeo$freqDif[chr4_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[1,4] <- pval

DiavAppleGeo_cors[1,4] <- corEst$estimate

#Chromosome 4 High LD

SD_mat <- read.table("./data/SD_genosChr4_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDH_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr4_LDH_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr4_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_H_ind_p], AppleGeo$freqDif[chr4_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[2,4] <- pval

DiavAppleGeo_cors[2,4] <- corEst$estimate

#Chromosome 4 int LD

SD_mat <- read.table("./data/SD_genosChr4_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDM_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr4_LDM_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr4_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_M_ind_p], AppleGeo$freqDif[chr4_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[3,4] <- pval

DiavAppleGeo_cors[3,4] <- corEst$estimate


#Chromosome 4 low LD

SD_mat <- read.table("./data/SD_genosChr4_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDL_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr4_LDL_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr4_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_L_ind_p], AppleGeo$freqDif[chr4_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[4,4] <- pval

DiavAppleGeo_cors[4,4] <- corEst$estimate

#Chromosome 5

SD_mat <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr5_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr5_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_ind_p], AppleGeo$freqDif[chr5_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[1,5] <- pval

DiavAppleGeo_cors[1,5] <- corEst$estimate

#Chromosome 5 High LD

SD_mat <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr5_LDH_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr5_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_H_ind_p], AppleGeo$freqDif[chr5_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[2,5] <- pval

DiavAppleGeo_cors[2,5] <- corEst$estimate

#Chromosome 5 int LD

SD_mat <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr5_LDM_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr5_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_M_ind_p], AppleGeo$freqDif[chr5_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[3,5] <- pval

DiavAppleGeo_cors[3,5] <- corEst$estimate


#Chromosome 5 low LD

SD_mat <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genosChr5_LDL_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genosChr5_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_L_ind_p], AppleGeo$freqDif[chr5_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[4,5] <- pval

DiavAppleGeo_cors[4,5] <- corEst$estimate


#All Mapped loci 
SD_mat <- read.table("./data/SD_genosAllChr_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosAllChr_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genos_AllChr_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genos_AllChr_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_chr_ind_p], AppleGeo$freqDif[all_chr_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[1,6] <- pval

DiavAppleGeo_cors[1,6] <- corEst$estimate

#All high LD 
SD_mat <- read.table("./data/SD_genos_AllLDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDH_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genos_AllLDH_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genos_AllLDH_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDH_p], AppleGeo$freqDif[all_LDH_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[2,6] <- pval

DiavAppleGeo_cors[2,6] <- corEst$estimate

#All int. LD 
SD_mat <- read.table("./data/SD_genos_AllLDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDM_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genos_AllLDM_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genos_AllLDM_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDM_p], AppleGeo$freqDif[all_LDM_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[3,6] <- pval

DiavAppleGeo_cors[3,6] <- corEst$estimate


#All low LD
SD_mat <- read.table("./data/SD_genos_AllLDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDL_polarized.txt", header=T)
apple_7 <- read.table("./data/apple7_genos_AllLDL_polarized.txt", header=T)
UrbApple <- read.table("./data/UrbApple_genos_AllLDL_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDL_p], AppleGeo$freqDif[all_LDL_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_7, ecl2 = UrbApple, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavAppleGeo_pvals[4,6] <- pval

DiavAppleGeo_cors[4,6] <- corEst$estimate

# Write the estimate and pvalue files
write.table(DiavAppleGeo_cors, file = "./results/DiavAppleGeo_corEst.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(DiavAppleGeo_pvals, file = "./results/DiavAppleGeo_perPvals.txt", row.names = T, col.names = T, quote = F, sep = "\t")
#######################################
# Figures
#######################################

#Figure 5 and supplementary figure 3

#Chromosome 1 

dat <- as.data.frame(cbind(AppleGeo$freqDif, DiaInt$freqDif))
lowld <- as.data.frame(cbind(AppleGeo$freqDif[chr1_L_ind_p], DiaInt$freqDif[chr1_L_ind_p]))
midld <- as.data.frame(cbind(AppleGeo$freqDif[chr1_M_ind_p], DiaInt$freqDif[chr1_M_ind_p]))
highld <- as.data.frame(cbind(AppleGeo$freqDif[chr1_H_ind_p], DiaInt$freqDif[chr1_H_ind_p]))
colnames(dat) <- c("AppleGeo", "Dia")
colnames(lowld) <- c("AppleGeo", "Dia")
colnames(midld) <- c("AppleGeo", "Dia")
colnames(highld) <- c("AppleGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavAppleGeo_chr1.tiff")
p <- ggplot(dat, aes(x=AppleGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=AppleGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=AppleGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=AppleGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Apple Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 2 
lowld <- as.data.frame(cbind(AppleGeo$freqDif[chr2_L_ind_p], DiaInt$freqDif[chr2_L_ind_p]))
midld <- as.data.frame(cbind(AppleGeo$freqDif[chr2_M_ind_p], DiaInt$freqDif[chr2_M_ind_p]))
highld <- as.data.frame(cbind(AppleGeo$freqDif[chr2_H_ind_p], DiaInt$freqDif[chr2_H_ind_p]))
colnames(dat) <- c("AppleGeo", "Dia")
colnames(lowld) <- c("AppleGeo", "Dia")
colnames(midld) <- c("AppleGeo", "Dia")
colnames(highld) <- c("AppleGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavAppleGeo_chr2.tiff")
p <- ggplot(dat, aes(x=AppleGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=AppleGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=AppleGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=AppleGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Apple Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 3 
lowld <- as.data.frame(cbind(AppleGeo$freqDif[chr3_L_ind_p], DiaInt$freqDif[chr3_L_ind_p]))
midld <- as.data.frame(cbind(AppleGeo$freqDif[chr3_M_ind_p], DiaInt$freqDif[chr3_M_ind_p]))
highld <- as.data.frame(cbind(AppleGeo$freqDif[chr3_H_ind_p], DiaInt$freqDif[chr3_H_ind_p]))
colnames(dat) <- c("AppleGeo", "Dia")
colnames(lowld) <- c("AppleGeo", "Dia")
colnames(midld) <- c("AppleGeo", "Dia")
colnames(highld) <- c("AppleGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavAppleGeo_chr3.tiff")
p <- ggplot(dat, aes(x=AppleGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=AppleGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=AppleGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=AppleGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Apple Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 4
lowld <- as.data.frame(cbind(AppleGeo$freqDif[chr4_L_ind_p], DiaInt$freqDif[chr4_L_ind_p]))
midld <- as.data.frame(cbind(AppleGeo$freqDif[chr4_M_ind_p], DiaInt$freqDif[chr4_M_ind_p]))
highld <- as.data.frame(cbind(AppleGeo$freqDif[chr4_H_ind_p], DiaInt$freqDif[chr4_H_ind_p]))
colnames(dat) <- c("AppleGeo", "Dia")
colnames(lowld) <- c("AppleGeo", "Dia")
colnames(midld) <- c("AppleGeo", "Dia")
colnames(highld) <- c("AppleGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavAppleGeo_chr4.tiff")
p <- ggplot(dat, aes(x=AppleGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=AppleGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=AppleGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=AppleGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Apple Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 5
lowld <- as.data.frame(cbind(AppleGeo$freqDif[chr5_L_ind_p], DiaInt$freqDif[chr5_L_ind_p]))
midld <- as.data.frame(cbind(AppleGeo$freqDif[chr5_M_ind_p], DiaInt$freqDif[chr5_M_ind_p]))
highld <- as.data.frame(cbind(AppleGeo$freqDif[chr5_H_ind_p], DiaInt$freqDif[chr5_H_ind_p]))
colnames(dat) <- c("AppleGeo", "Dia")
colnames(lowld) <- c("AppleGeo", "Dia")
colnames(midld) <- c("AppleGeo", "Dia")
colnames(highld) <- c("AppleGeo", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavAppleGeo_chr5.tiff")
p <- ggplot(dat, aes(x=AppleGeo, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=AppleGeo, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=AppleGeo, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=AppleGeo, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Apple Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

