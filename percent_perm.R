#Percent of significant SNPs permutation test. 

# This analysis tests whether the number of significant SNPs in a given chromosome - LD group is more than 
# would be expected by chance. 

#To perform this test, A large matrix of permuted allele frequency differences is created. With columns 
# being loci and rows different iterations of the permutation. After the permuted allele frequency differences 
# are calculated, they are converted to what percentile of the distribution of permutations they fall on. Then
# for each locus in a given replication, the number of significant SNPs (< 0.05 cut off) is summed and divided by 
# the number of loci in the row (i.e. the chromosme-LD group under study). These percents are calculated for every 
# permutation to create a null distribution of the percent of significant SNPs that would be expected on given group. 

setwd("/media/raglandlab/ExtraDrive4/Clines_3/")

#Load in R functions
source("./src/clines3_functions.R")

library(iterators)
library(doParallel)
library(foreach)

#Load in the mapping information: 
load("./data/Mapped_RAD_loci.Rdata")


#allele freq difference 
allele_freq <- read.table("./results/SDvCDND_freqDifs.txt", header=T)

################
# Chromosome 1 #
################
#LDL
#Load in data
SD_mat <- read.table("./data/SD_genosChr1_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDL_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)
#Convert to pvals 
cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)

#Get number of significant loci for each iteration
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr1_L_ind_p] < 0.05)/length(allele_freq$p_val[chr1_L_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- c(emp_perc, emp_perc_pval)
names(sum_r) <- c("emp_perc", "p_val") 
  
#LDM
#Load in data
SD_mat <- read.table("./data/SD_genosChr1_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDM_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr1_M_ind_p] < 0.05)/length(allele_freq$p_val[chr1_M_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


#LDH
#Load in data
SD_mat <- read.table("./data/SD_genosChr1_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr1_LDH_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr1_H_ind_p] < 0.05)/length(allele_freq$p_val[chr1_H_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


################
# Chromosome 2 #
################

#LDL
#Load in data
SD_mat <- read.table("./data/SD_genosChr2_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDL_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)
#Convert to pvals
cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)

#Get number of significant loci for each iteration
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr2_L_ind_p] < 0.05)/length(allele_freq$p_val[chr2_L_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- c(emp_perc, emp_perc_pval)
names(sum_r) <- c("emp_perc", "p_val")

#LDM
#Load in data
SD_mat <- read.table("./data/SD_genosChr2_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDM_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr2_M_ind_p] < 0.05)/length(allele_freq$p_val[chr2_M_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


#LDH
#Load in data
SD_mat <- read.table("./data/SD_genosChr2_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr2_LDH_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr2_H_ind_p] < 0.05)/length(allele_freq$p_val[chr2_H_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


################
# Chromosome 3 #
################

#LDL
#Load in data
SD_mat <- read.table("./data/SD_genosChr3_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDL_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)
#Convert to pvals
cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)

#Get number of significant loci for each iteration
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr3_L_ind_p] < 0.05)/length(allele_freq$p_val[chr3_L_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- c(emp_perc, emp_perc_pval)
names(sum_r) <- c("emp_perc", "p_val")

#LDM
#Load in data
SD_mat <- read.table("./data/SD_genosChr3_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDM_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr3_M_ind_p] < 0.05)/length(allele_freq$p_val[chr3_M_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


#LDH
#Load in data
SD_mat <- read.table("./data/SD_genosChr3_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr3_LDH_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr3_H_ind_p] < 0.05)/length(allele_freq$p_val[chr3_H_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


################
# Chromosome 4 #
################

#LDL
#Load in data
SD_mat <- read.table("./data/SD_genosChr4_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDL_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)
#Convert to pvals
cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)

#Get number of significant loci for each iteration
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr4_L_ind_p] < 0.05)/length(allele_freq$p_val[chr4_L_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- c(emp_perc, emp_perc_pval)
names(sum_r) <- c("emp_perc", "p_val")

#LDM
#Load in data
SD_mat <- read.table("./data/SD_genosChr4_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDM_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr4_M_ind_p] < 0.05)/length(allele_freq$p_val[chr4_M_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


#LDH
#Load in data
SD_mat <- read.table("./data/SD_genosChr4_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr4_LDH_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr4_H_ind_p] < 0.05)/length(allele_freq$p_val[chr4_H_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


################
# Chromosome 5 #
################

#LDL
#Load in data
SD_mat <- read.table("./data/SD_genosChr5_LDL_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDL_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)
#Convert to pvals
cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)

#Get number of significant loci for each iteration
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr5_L_ind_p] < 0.05)/length(allele_freq$p_val[chr5_L_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- c(emp_perc, emp_perc_pval)
names(sum_r) <- c("emp_perc", "p_val")

#LDM
#Load in data
SD_mat <- read.table("./data/SD_genosChr5_LDM_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDM_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr5_M_ind_p] < 0.05)/length(allele_freq$p_val[chr5_M_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))


#LDH
#Load in data
SD_mat <- read.table("./data/SD_genosChr5_LDH_polarized.txt", header=T)
CDND_mat <- read.table("./data/CDND_genosChr5_LDH_polarized.txt", header=T)

#register 10 cores
cl <- makeCluster(10)
registerDoParallel(cl)
nloci <- nrow(SD_mat)
#Permutation test
perm_mat<-foreach(i=1:nloci, .combine=cbind) %dopar% {
  permTest_perc(dia1 = unlist(SD_mat[i,]), dia2 = unlist(CDND_mat[i,]), nit = 10000)
}
stopCluster(cl)

cl <- makeCluster(10)
registerDoParallel(cl)
perm_pval_mat<-foreach(i=1:nrow(perm_mat), .combine='rbind') %:%
  foreach(j=1:ncol(perm_mat), .combine='c') %dopar% {
    permPval(perm_mat[i,j], perm_mat[,j])
  }
stopCluster(cl)
percs <- vector()
for(i in 1:nrow(perm_pval_mat)){
  percs[i] <- sum(perm_pval_mat[i,] < 0.05)/length(perm_pval_mat[i,])
}
emp_perc <- sum(allele_freq$p_val[chr5_H_ind_p] < 0.05)/length(allele_freq$p_val[chr5_H_ind_p])
emp_perc_pval <- permPval(emp_perc, percs)
sum_r <- rbind(sum_r,c(emp_perc, emp_perc_pval))

save.image(file = "percents.RData")

write.table(sum_r, file="./results/percent_permutation.txt", quote=F, row.names=F, col.names=T, sep="\t")

#Ok, now we'll graph this data. 
# chrs <- c(rep("1",3), rep("2", 3), rep("3",3), rep("4", 3), rep("5",3))
# LD <-c(rep(c("high", "int", "low"), 5))
# percent_df <- as.data.frame(cbind(percents, chrs, LD))
# colnames(percent_df) <- c("percents", "chromosomes", "LD")
# percent_df$percents <- as.numeric(as.character(percent_df$percents))
# #val <- c("", "", "", "****", "****", "***","****", "****", "", "", "", "", "****", "****", "*" )
# val <- c("", "", "", "*", "*", "*","*", "*", "", "", "", "", "*", "*", "*" )
# library(ggplot2)
# library(RColorBrewer)
# library(viridis)
# pdf("CDND-SD_percents.pdf")
# mypalette <- viridis_pal(option = "B")(1000)
# p<- ggplot(data=percent_df, aes(x=chromosomes, y=percents, fill=LD))+
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   theme_classic() + scale_y_continuous(limits = c(0,1))+ xlab("chromosomes") +
#   ylab("percent significant snps")
# p + scale_fill_manual(values = c(mypalette[250], mypalette[500], mypalette[750])) + 
#   geom_text(aes(label = val), vjust = -0.5, position = position_dodge(0.9), size=16)
# #theme(axis.line = element_line(size = 2), axis.text = element_text(size = 24.0))
# dev.off()


