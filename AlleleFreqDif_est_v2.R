#Allele frequency difference estimates for each of the datasets in this study. 

#These scripts take advantage of Parrellization and are designed to be run on a comptuer with multiple 
#processors or a high performance computing cluster. 

#They may be performed on smaller computers, such as a laptop, but the code will need to be adjusted accordingly.

#1. Diapause class (This study)
#2. Selection experiments (Egan et al. 2015 [haw host race], Doellman et al. 2019 [apple host race])
#3. Diapause termination GWAS (Ragland et al. 2017)
#4. Geographic variation (Doellman et al. 2018)

setwd("/media/raglandlab/ExtraDrive4/Clines_3/")

#Load in R functions
source("./src/clines3_functions.R")

library(iterators)
library(doParallel)
library(foreach)

#Load in data
SD_mat <- read.table("./data/SD_genosAll_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosAll_polarized.txt", header=T)

#Load in the chromosome 5 data 
SD_mat_chr5 <- read.table("./data/SD_genosChr5_polarized.txt", header=T)
CDND_mat_chr5 <- read.table("./data/CDND_genosChr5_polarized.txt", header = T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Permutation test
pval_shallowVsCDND<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(SD_mat[i,]),unlist(CDND_mat[i,]),nit=10000)
}

#Empirical Frequeny difference estimate
pointEst_shallowVsCDND<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(SD_mat[i,]),unlist(CDND_mat[i,]))
}

#now for the chromosome 5 data 
nloci=nrow(SD_mat_chr5)
pval_shallowVsCDND_chr5 <-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(SD_mat_chr5[i,]),unlist(CDND_mat_chr5[i,]),nit=10000)
}

pointEst_shallowVsCDND_chr5<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(SD_mat_chr5[i,]),unlist(CDND_mat_chr5[i,]))
}

stopCluster(cl)

#bind up the frequency differenes etimates and p values and print the file. 
CDNDvSD <- cbind(pointEst_shallowVsCDND, pval_shallowVsCDND)
rownames(CDNDvSD) <- rownames(SD_mat)
colnames(CDNDvSD) <- c("freqDif", "p_val")

#bind up the chromosome 5 data 
CDNDvSD_chr5 <- cbind(pointEst_shallowVsCDND_chr5, pval_shallowVsCDND_chr5)
rownames(CDNDvSD_chr5) <- rownames(SD_mat_chr5)
colnames(CDNDvSD_chr5) <- c("freqDif", "p_val")

#Replace chromosome 5 loci in complete set (CDNDvSD) with the female only calculate of chromosome 5 
ind <- rownames(CDNDvSD) %in% rownames(CDNDvSD_chr5)
CDNDvSD[ind,] <- CDNDvSD_chr5

write.table(CDNDvSD, file="./results/SDvCDND_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#Now, we'll do it pairiwse for SD - CD, SD - ND, and CD-ND 

CD_mat <- read.table("./data/CD_genosAll_polarized.txt", header=T)
ND_mat <- read.table("./data/ND_genosAll_polarized.txt", header=T)

#Chromosome 5 
CD_mat_chr5 <- read.table("./data/CD_genosChr5_polarized.txt", header=T)
ND_mat_chr5 <- read.table("./data/ND_genosChr5_polarized.txt", header=T)
#SD vs. CD

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_SDvCD<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(SD_mat[i,]),unlist(CD_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_SDvCD<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(SD_mat[i,]),unlist(CD_mat[i,]))
}

#now for the chromosome 5 data 
nloci=nrow(SD_mat_chr5)
pval_SDvCD_chr5 <-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(SD_mat_chr5[i,]),unlist(CD_mat_chr5[i,]),nit=10000)
}

pointEst_SDvCD_chr5<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(SD_mat_chr5[i,]),unlist(CD_mat_chr5[i,]))
}

stopCluster(cl)

SDvCD <- cbind(pointEst_SDvCD, pval_SDvCD)
rownames(SDvCD) <- rownames(SD_mat)
colnames(SDvCD) <- c("freqDif", "p_val")

#bind up the chromosome 5 data 
SDvCD_chr5 <- cbind(pointEst_SDvCD_chr5, pval_SDvCD_chr5)
rownames(SDvCD_chr5) <- rownames(SD_mat_chr5)
colnames(SDvCD_chr5) <- c("freqDif", "p_val")

#Replace chromosome 5 loci in complete set with the female only calculate of chromosome 5 
ind <- rownames(SDvCD) %in% rownames(SDvCD_chr5)
SDvCD[ind,] <- SDvCD_chr5

write.table(CDNDvSD, file="./results/SDvCD_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#SD v ND 

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_SDvND<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(SD_mat[i,]),unlist(ND_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_SDvND<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(SD_mat[i,]),unlist(ND_mat[i,]))
}

nloci=nrow(SD_mat_chr5)
pval_SDvND_chr5 <-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(SD_mat_chr5[i,]),unlist(ND_mat_chr5[i,]),nit=10000)
}

pointEst_SDvND_chr5<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(SD_mat_chr5[i,]),unlist(ND_mat_chr5[i,]))
}

stopCluster(cl)

SDvND <- cbind(pointEst_SDvND, pval_SDvND)
rownames(SDvND) <- rownames(SD_mat)
colnames(SDvND) <- c("freqDif", "p_val")

#bind up the chromosome 5 data 
SDvND_chr5 <- cbind(pointEst_SDvND_chr5, pval_SDvND_chr5)
rownames(SDvND_chr5) <- rownames(SD_mat_chr5)
colnames(SDvND_chr5) <- c("freqDif", "p_val")

#Replace chromosome 5 loci in complete set with the female only calculate of chromosome 5 
ind <- rownames(SDvCD) %in% rownames(SDvND_chr5)
SDvCD[ind,] <- SDvND_chr5


write.table(SDvND, file="./results/SDvND_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#CD vs. ND

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_CDvND<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(CD_mat[i,]),unlist(ND_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_CDvND<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(CD_mat[i,]),unlist(ND_mat[i,]))
}

nloci=nrow(CD_mat_chr5)
pval_CDvND_chr5 <-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(CD_mat_chr5[i,]),unlist(ND_mat_chr5[i,]),nit=10000)
}

pointEst_CDvND_chr5<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(CD_mat_chr5[i,]),unlist(ND_mat_chr5[i,]))
}

stopCluster(cl)

CDvND <- cbind(pointEst_CDvND, pval_CDvND)
rownames(CDvND) <- rownames(CD_mat)
colnames(CDvND) <- c("freqDif", "p_val")

CDvND_chr5 <- cbind(pointEst_CDvND_chr5, pval_CDvND_chr5)
rownames(CDvND_chr5) <- rownames(CD_mat_chr5)
colnames(CDvND_chr5) <- c("freqDif", "p_val")

#Replace chromosome 5 loci in complete set with the female only calculate of chromosome 5 
ind <- rownames(CDvND) %in% rownames(CDvND_chr5)
CDvND[ind,] <- CDvND_chr5


write.table(CDvND, file="./results/CDvND_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")


####################################
# 2. AFD for Selection experiments 
###################################

## Apple - haw for Grant comparison 

apple_7_mat <- read.table("./data/apple7_genosAll_polarized.txt", header=T)
haw_7_mat <- read.table("./data/haw7_genosAll_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_apple_7vhaw_7<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(apple_7_mat[i,]),unlist(haw_7_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_apple_7vhaw_7<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(apple_7_mat[i,]),unlist(haw_7_mat[i,]))
}
stopCluster(cl)

apple_7vhaw_7 <- cbind(pointEst_apple_7vhaw_7, pval_apple_7vhaw_7)
rownames(apple_7vhaw_7) <- rownames(apple_7_mat)
colnames(apple_7vhaw_7) <- c("freqDif", "p_val")
write.table(apple_7vhaw_7, file="./results/apple7haw7_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#Apple selection experiment

apple_7_mat <- read.table("./data/apple7_genosAll_polarized.txt", header=T)
apple_32_mat <- read.table("./data/apple32_genosAll_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_apple_7vapple_32<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(apple_7_mat[i,]),unlist(apple_32_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_apple_7vapple_32<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(apple_7_mat[i,]),unlist(apple_32_mat[i,]))
}
stopCluster(cl)

apple_7vapple_32 <- cbind(pointEst_apple_7vapple_32, pval_apple_7vapple_32)
rownames(apple_7vapple_32) <- rownames(apple_7_mat)
colnames(apple_7vapple_32) <- c("freqDif", "p_val")
write.table(apple_7vapple_32, file="./results/appleSel_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#Haw selection experiment

haw_7_mat <- read.table("./data/haw7_genosAll_polarized.txt", header=T)
haw_32_mat <- read.table("./data/haw32_genosAll_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_haw_7vhaw_32<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(haw_7_mat[i,]),unlist(haw_32_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_haw_7vhaw_32<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(haw_7_mat[i,]),unlist(haw_32_mat[i,]))
}
stopCluster(cl)

haw_7vhaw_32 <- cbind(pointEst_haw_7vhaw_32, pval_haw_7vhaw_32)
rownames(haw_7vhaw_32) <- rownames(haw_7_mat)
colnames(haw_7vhaw_32) <- c("freqDif", "p_val")
write.table(haw_7vhaw_32, file="./results/hawSel_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

####################################
# 3. Diapaue termination
###################################

#Apple race 

apple_early_mat <- read.table("./data/appleEarly_genosAll_polarized.txt", header=T)
apple_late_mat <- read.table("./data/appleLate_genosAll_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_apple_earlyvapple_late<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(apple_early_mat[i,]),unlist(apple_late_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEstimate_apple_earlyvapple_late<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(apple_early_mat[i,]),unlist(apple_late_mat[i,]))
}
stopCluster(cl)

apple_earlyvapple_late <- cbind(pointEstimate_apple_earlyvapple_late, pval_apple_earlyvapple_late)
rownames(apple_earlyvapple_late) <- rownames(apple_early_mat)
colnames(apple_earlyvapple_late) <- c("freqDif", "p_val")
write.table(apple_earlyvapple_late, file="./results/apple_EarlyvLate_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")


#Haw race 

haw_early_mat <- read.table("./data/hawEarly_genosAll_polarized.txt", header=T)
haw_late_mat <- read.table("./data/hawLate_genosAll_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_haw_earlyvhaw_late<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(haw_early_mat[i,]),unlist(haw_late_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEstimate_haw_earlyvhaw_late<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(haw_early_mat[i,]),unlist(haw_late_mat[i,]))
}
stopCluster(cl)

haw_earlyvhaw_late <- cbind(pointEstimate_haw_earlyvhaw_late, pval_haw_earlyvhaw_late)
rownames(haw_earlyvhaw_late) <- rownames(haw_early_mat)
colnames(haw_earlyvhaw_late) <- c("freqDif", "p_val")
write.table(haw_earlyvhaw_late, file="./results/haw_EarlyvLate_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#Average early - late 

#Data has already been loaded in above. 

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_avg_early_late<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest_2(unlist(apple_early_mat[i,]),unlist(apple_late_mat[i,]), unlist(haw_early_mat[i,]), unlist(haw_late_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_avg_early_late<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif_2(f1 = unlist(apple_early_mat[i,]), f2 = unlist(apple_late_mat[i,]), f3 = unlist(haw_early_mat[i,]), f4 = unlist(haw_late_mat[i,]))
}
stopCluster(cl)

avg_early_late <- cbind(pointEst_avg_early_late, pval_avg_early_late)
rownames(avg_early_late) <- rownames(haw_early_mat)
colnames(avg_early_late) <- c("freqDif", "p_val")
write.table(avg_early_late, file="./results/avg_early_late_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

####################################
# 4. Geographic variation
###################################

#This estimate allele frequency differences between our northern most (Grant, MI) and southern most (Urbana, IL) sites.

#Apple race

apple_Gra_mat <- read.table("./data/apple7_genosAll_polarized.txt", header=T)
apple_Urb_mat <- read.table("./data/UrbApple_genosAllChr_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_apple_Gravapple_Urb<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(apple_Gra_mat[i,]),unlist(apple_Urb_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_apple_Gravapple_Urb<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(apple_Gra_mat[i,]),unlist(apple_Urb_mat[i,]))
}
stopCluster(cl)

apple_Gravapple_Urb <- cbind(pointEst_apple_Gravapple_Urb, pval_apple_Gravapple_Urb)
rownames(apple_Gravapple_Urb) <- rownames(apple_Gra_mat)
colnames(apple_Gravapple_Urb) <- c("freqDif", "p_val")
write.table(apple_Gravapple_Urb, file="./results/Apple_geo_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")


#Haw race

haw_Gra_mat <- read.table("./data/haw7_genosAll_polarized.txt", header=T)
haw_Urb_mat <- read.table("./data/UrbHaw_genosAllChr_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
#Produce permutation matrix
pval_haw_Gravhaw_Urb<-foreach(i=1:nloci, .combine=c) %dopar% {
  permTest(unlist(haw_Gra_mat[i,]),unlist(haw_Urb_mat[i,]),nit=10000)
}
#Frequeny difference estimate
pointEst_haw_Gravhaw_Urb<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqDif(unlist(haw_Gra_mat[i,]),unlist(haw_Urb_mat[i,]))
}
stopCluster(cl)

haw_Gravhaw_Urb <- cbind(pointEst_haw_Gravhaw_Urb, pval_haw_Gravhaw_Urb)
rownames(haw_Gravhaw_Urb) <- rownames(haw_Gra_mat)
colnames(haw_Gravhaw_Urb) <- c("freqDif", "p_val")
write.table(haw_Gravhaw_Urb, file="./results/Haw_geo_freqDifs.txt", row.names = T, col.names = T, quote = F, sep = "\t")

