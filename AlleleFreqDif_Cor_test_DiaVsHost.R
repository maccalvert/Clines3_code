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

##################################################
# 1. Diapause intensity vs. Host race (Grant, MI)#
##################################################

#Load in the empirical allele frequency differences for each 
DiaInt <- read.table("./results/SDvCDND_freqDifs.txt", header=T)
HostDif <- read.table("./results/apple7haw7_freqDifs.txt", header=T)

#Make empty matrices to store the p_values and corelation estimates 
DiavHost_pvals <- matrix(nrow=4, ncol =6)
rownames(DiavHost_pvals) <- c("All_Snps", "High_LD", "Int. LD", "Low LD")
colnames(DiavHost_pvals) <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_all")
DiavHost_cors <- matrix(nrow=4, ncol =6)
rownames(DiavHost_cors) <- c("All_Snps", "High_LD", "Int. LD", "Low LD")
colnames(DiavHost_cors) <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_all")

#Chromosome 1 

SD_mat <- read.table("./data/SD_genosChr1_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr1_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr1_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_ind_p], HostDif$freqDif[chr1_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[1,1] <- pval

DiavHost_cors[1,1] <- corEst$estimate

#Chromosome 1 High LD

SD_mat <- read.table("./data/SD_genosChr1_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDH_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr1_LDH_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr1_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_H_ind_p], HostDif$freqDif[chr1_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[2,1] <- pval

DiavHost_cors[2,1] <- corEst$estimate

#Chromosome 1 Int. LD 

SD_mat <- read.table("./data/SD_genosChr1_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDM_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr1_LDM_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr1_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_M_ind_p], HostDif$freqDif[chr1_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[3,1] <- pval

DiavHost_cors[3,1] <- corEst$estimate

#Chromosome 1 Low LD 

SD_mat <- read.table("./data/SD_genosChr1_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDL_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr1_LDL_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr1_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr1_L_ind_p], HostDif$freqDif[chr1_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[4,1] <- pval

DiavHost_cors[4,1] <- corEst$estimate

#Chromosome 2

SD_mat <- read.table("./data/SD_genosChr2_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr2_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr2_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_ind_p], HostDif$freqDif[chr2_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[1,2] <- pval

DiavHost_cors[1,2] <- corEst$estimate

#Chromosome 2 High LD

SD_mat <- read.table("./data/SD_genosChr2_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDH_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr2_LDH_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr2_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_H_ind_p], HostDif$freqDif[chr2_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[2,2] <- pval

DiavHost_cors[2,2] <- corEst$estimate

#Chromosome 2 Int. LD 

SD_mat <- read.table("./data/SD_genosChr2_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDM_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr2_LDM_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr2_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_M_ind_p], HostDif$freqDif[chr2_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[3,2] <- pval

DiavHost_cors[3,2] <- corEst$estimate

#Chromosome 2 low LD 

SD_mat <- read.table("./data/SD_genosChr2_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDL_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr2_LDL_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr2_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr2_L_ind_p], HostDif$freqDif[chr2_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[4,2] <- pval

DiavHost_cors[4,2] <- corEst$estimate

#Chromosome 3

SD_mat <- read.table("./data/SD_genosChr3_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr3_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr3_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_ind_p], HostDif$freqDif[chr3_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[1,3] <- pval

DiavHost_cors[1,3] <- corEst$estimate

#Chromosome 3 High LD

SD_mat <- read.table("./data/SD_genosChr3_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDH_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr3_LDH_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr3_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_H_ind_p], HostDif$freqDif[chr3_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[2,3] <- pval

DiavHost_cors[2,3] <- corEst$estimate

#Chromosome 3 Int. LD 

SD_mat <- read.table("./data/SD_genosChr3_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDM_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr3_LDM_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr3_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_M_ind_p], HostDif$freqDif[chr3_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[3,3] <- pval

DiavHost_cors[3,3] <- corEst$estimate

#Chromosome 3 Low LD 

SD_mat <- read.table("./data/SD_genosChr3_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDL_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr3_LDL_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr3_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr3_L_ind_p], HostDif$freqDif[chr3_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[4,3] <- pval

DiavHost_cors[4,3] <- corEst$estimate

#Chromosome 4

SD_mat <- read.table("./data/SD_genosChr4_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr4_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr4_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_ind_p], HostDif$freqDif[chr4_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[1,4] <- pval

DiavHost_cors[1,4] <- corEst$estimate

#Chromosome 4 High LD

SD_mat <- read.table("./data/SD_genosChr4_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDH_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr4_LDH_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr4_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_H_ind_p], HostDif$freqDif[chr4_H_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[2,4] <- pval

DiavHost_cors[2,4] <- corEst$estimate

#Chromosome 4 Int. LD 

SD_mat <- read.table("./data/SD_genosChr4_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDM_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr4_LDM_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr4_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_M_ind_p], HostDif$freqDif[chr4_M_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[3,4] <- pval

DiavHost_cors[3,4] <- corEst$estimate

#Chromosome 4 Low LD 

SD_mat <- read.table("./data/SD_genosChr4_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDL_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr4_LDL_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr4_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr4_L_ind_p], HostDif$freqDif[chr4_L_ind_p], method = "spearman")

cl <- makeCluster(10)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[4,4] <- pval

DiavHost_cors[4,4] <- corEst$estimate

#Chromosome 5

SD_mat <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr5_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr5_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_ind_p], HostDif$freqDif[chr5_ind_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors<-foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[1,5] <- pval

DiavHost_cors[1,5] <- corEst$estimate

#Chromosome 5 High LD

SD_mat <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr5_LDH_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr5_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_H_ind_p], HostDif$freqDif[chr5_H_ind_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[2,5] <- pval

DiavHost_cors[2,5] <- corEst$estimate

#Chromosome 5 Int. LD 

SD_mat <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr5_LDM_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr5_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_M_ind_p], HostDif$freqDif[chr5_M_ind_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[3,5] <- pval

DiavHost_cors[3,5] <- corEst$estimate


#Chromosome 5 Low LD 
SD_mat <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genosChr5_LDL_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genosChr5_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[chr5_L_ind_p], HostDif$freqDif[chr5_L_ind_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine=c) %dopar% {
  permCor(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[4,5] <- pval

DiavHost_cors[4,5] <- corEst$estimate

#All mapped chromosomes 
SD_mat <- read.table("./data/SD_genosAllChr_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosAllChr_polarized.txt", header=T)
apple_mat <- read.table("./data/apple7_genos_AllChr_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genos_AllChr_polarized.txt", header=T)

#Chr5 loci 
SD_mat_chr5 <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_chr_ind_p], HostDif$freqDif[all_chr_ind_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[1,6] <- pval

DiavHost_cors[1,6] <- corEst$estimate

#All High LD 
SD_mat <- read.table("./data/SD_genos_AllLDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDH_polarized.txt", header=T)
apple_mat <- read.table("./data/haw7_genos_AllLDH_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genos_AllLDH_polarized.txt", header=T)

#Chromosome 5 
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDH_p], HostDif$freqDif[all_LDH_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[2,6] <- pval

DiavHost_cors[2,6] <- corEst$estimate

#All Int. LD 
SD_mat <- read.table("./data/SD_genos_AllLDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDM_polarized.txt", header=T)
apple_mat <- read.table("./data/haw7_genos_AllLDM_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genos_AllLDM_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDM_p], HostDif$freqDif[all_LDM_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[3,6] <- pval

DiavHost_cors[3,6] <- corEst$estimate

#All Low LD
SD_mat <- read.table("./data/SD_genos_AllLDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_AllLDL_polarized.txt", header=T)
apple_mat <- read.table("./data/haw7_genos_AllLDL_polarized.txt", header=T)
haw_mat <- read.table("./data/haw7_genos_AllLDL_polarized.txt", header=T)

#Chromosome 5
SD_mat_chr5 <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)

corEst <- cor.test(DiaInt$freqDif[all_LDL_p], HostDif$freqDif[all_LDL_p], method = "spearman")

cl <- makeCluster(8)
registerDoParallel(cl)
permuted_Cors <- foreach(i=1:10, .combine = c) %dopar% {
  permCorForLDs(dia1 = SD_mat, dia2 = CDND_mat, ecl1 = apple_mat, ecl2 = haw_mat, chr5_1 = SD_mat_chr5, chr5_2 = CDND_mat_chr5, nit=((i-(i-1))*1000))
}
stopCluster(cl)

pval <- permPval(corEst$estimate, permuted_Cors)

DiavHost_pvals[4,6] <- pval

DiavHost_cors[4,6] <- corEst$estimate

# Write the estimate and pvalue files
write.table(DiavHost_cors, file = "./results/DiavHost_corEst.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(DiavHost_pvals, file = "./results/DiavHost_perPvals.txt", row.names = T, col.names = T, quote = F, sep = "\t")
#######################################
# Figures
#######################################

#Figure 5 and supplementary figure 3

#Chromosome 1 

dat <- as.data.frame(cbind(HostDif$freqDif, DiaInt$freqDif))
lowld <- as.data.frame(cbind(HostDif$freqDif[chr1_L_ind_p], DiaInt$freqDif[chr1_L_ind_p]))
midld <- as.data.frame(cbind(HostDif$freqDif[chr1_M_ind_p], DiaInt$freqDif[chr1_M_ind_p]))
highld <- as.data.frame(cbind(HostDif$freqDif[chr1_H_ind_p], DiaInt$freqDif[chr1_H_ind_p]))
colnames(dat) <- c("HostRace", "Dia")
colnames(lowld) <- c("HostRace", "Dia")
colnames(midld) <- c("HostRace", "Dia")
colnames(highld) <- c("HostRace", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHost_chr1.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HostRace, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HostRace, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HostRace, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Countour plots
#pdf("./results/HostVint_cont_chr1.pdf")
tiff("./results/HostVint_cont_chr1.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(color="black", data = dat[chr1_ind_p,], aes(x=HostRace, y=Dia))+
  geom_density_2d(data = dat[chr1_ind_p,], aes(x=HostRace, y=Dia), binwidth=3)+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 2 
lowld <- as.data.frame(cbind(HostDif$freqDif[chr2_L_ind_p], DiaInt$freqDif[chr2_L_ind_p]))
midld <- as.data.frame(cbind(HostDif$freqDif[chr2_M_ind_p], DiaInt$freqDif[chr2_M_ind_p]))
highld <- as.data.frame(cbind(HostDif$freqDif[chr2_H_ind_p], DiaInt$freqDif[chr2_H_ind_p]))
colnames(dat) <- c("HostRace", "Dia")
colnames(lowld) <- c("HostRace", "Dia")
colnames(midld) <- c("HostRace", "Dia")
colnames(highld) <- c("HostRace", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHost_chr2.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HostRace, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HostRace, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HostRace, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Contour plots 
#pdf("./results/HostVint_cont_chr2.pdf")
tiff("./results/HostVint_cont_chr2.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(color="black", data = dat[chr2_ind_p,], aes(x=HostRace, y=Dia))+
  geom_density_2d(data = dat[chr2_ind_p,], aes(x=HostRace, y=Dia), binwidth=3)+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 3 
lowld <- as.data.frame(cbind(HostDif$freqDif[chr3_L_ind_p], DiaInt$freqDif[chr3_L_ind_p]))
midld <- as.data.frame(cbind(HostDif$freqDif[chr3_M_ind_p], DiaInt$freqDif[chr3_M_ind_p]))
highld <- as.data.frame(cbind(HostDif$freqDif[chr3_H_ind_p], DiaInt$freqDif[chr3_H_ind_p]))
colnames(dat) <- c("HostRace", "Dia")
colnames(lowld) <- c("HostRace", "Dia")
colnames(midld) <- c("HostRace", "Dia")
colnames(highld) <- c("HostRace", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHost_chr3.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HostRace, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HostRace, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HostRace, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Contour plots 
#pdf("./results/HostVint_cont_chr3.pdf")
tiff("./results/HostVint_cont_chr3.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(color="black", data = dat[chr3_ind_p,], aes(x=HostRace, y=Dia))+
  geom_density_2d(data = dat[chr3_ind_p,], aes(x=HostRace, y=Dia), binwidth=3)+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 4
lowld <- as.data.frame(cbind(HostDif$freqDif[chr4_L_ind_p], DiaInt$freqDif[chr4_L_ind_p]))
midld <- as.data.frame(cbind(HostDif$freqDif[chr4_M_ind_p], DiaInt$freqDif[chr4_M_ind_p]))
highld <- as.data.frame(cbind(HostDif$freqDif[chr4_H_ind_p], DiaInt$freqDif[chr4_H_ind_p]))
colnames(dat) <- c("HostRace", "Dia")
colnames(lowld) <- c("HostRace", "Dia")
colnames(midld) <- c("HostRace", "Dia")
colnames(highld) <- c("HostRace", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHost_chr4.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HostRace, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HostRace, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HostRace, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Contour plots 
#pdf("./results/HostVint_cont_chr4.pdf")
tiff("./results/HostVint_cont_chr4.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(color="black", data = dat[chr4_ind_p,], aes(x=HostRace, y=Dia))+
  geom_density_2d(data = dat[chr4_ind_p,], aes(x=HostRace, y=Dia), binwidth=3)+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Chromosome 5
lowld <- as.data.frame(cbind(HostDif$freqDif[chr5_L_ind_p], DiaInt$freqDif[chr5_L_ind_p]))
midld <- as.data.frame(cbind(HostDif$freqDif[chr5_M_ind_p], DiaInt$freqDif[chr5_M_ind_p]))
highld <- as.data.frame(cbind(HostDif$freqDif[chr5_H_ind_p], DiaInt$freqDif[chr5_H_ind_p]))
colnames(dat) <- c("HostRace", "Dia")
colnames(lowld) <- c("HostRace", "Dia")
colnames(midld) <- c("HostRace", "Dia")
colnames(highld) <- c("HostRace", "Dia")

mypallete <- viridis_pal(option = "B")(1000)
tiff("./results/DiavHost_chr5.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(data=lowld, aes(x=HostRace, y=Dia), color=mypallete[750])+
  geom_point(data=midld, aes(x=HostRace, y=Dia), color=mypallete[500])+
  geom_point(data=highld, aes(x=HostRace, y=Dia), color=mypallete[250])+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()

#Contour plots 
#pdf("./results/HostVint_cont_chr5.pdf")
tiff("./results/HostVint_cont_chr5.tiff")
p <- ggplot(dat, aes(x=HostRace, y=Dia))+
  geom_point(color="grey")+
  geom_point(color="black", data = dat[chr5_ind_p,], aes(x=HostRace, y=Dia))+
  geom_density_2d(data = dat[chr5_ind_p,], aes(x=HostRace, y=Dia), binwidth=3)+
  theme_classic()+
  theme(axis.text = element_text(size=24, color="black"))+
  ylab("SD - CD+ND Allele Freq. Dif.")+
  xlab("Apple - Haw Allele Freq. Dif.")+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = c(seq(-0.4,0.4,0.2)))+
  scale_x_continuous(limits = c(-0.05,0.4),breaks = c(seq(0,0.4,0.2)))
p
dev.off()
