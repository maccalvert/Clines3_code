#### These scripts analyze the results from the BSLMM , perform the polygenic score analyisis, and measure whether there is significant 
#overlap in the loci that are significantly differentiated between shallow diapause 
# many of the BSLMM scripts are adapted from those presented in Victor Soria-Carrasco's 
# Population Genomics workshop materials - https://github.com/visoca/popgenomworkshop-gwas_gemma

setwd("/media/raglandlab/ExtraDrive4/Clines_3")


#Overlap in loci significantly associated with diapause phenotypes. 
SDvCD <- read.table("./results/SDvCD_freqDifs.txt", header=T)
sum(SDvCD$p_val < 0.05)
#[1] 1139

SDvND <- read.table("./results/SDvND_freqDifs.txt", header=T)
sum(SDvND$p_val < 0.05)
#[1] 1177 

CDvND <- read.table("./results/CDvND_freqDifs.txt", header=T)
sum(CDvND$p_val < 0.05)
#[1] 394

#More loci than overlap by chance
# sig loci SD v CD    non sig loci SD v CD       
# a   b #sig loci SD v ND
# c   d #non sig loci SD v ND
a <- sum(SDvCD$p_val < 0.05 & SDvND$p_val < 0.05)
b <- sum(SDvCD$p_val > 0.05 & SDvND$p_val < 0.05)
c <- sum(SDvCD$p_val < 0.05 & SDvND$p_val > 0.05)
d <- sum(SDvCD$p_val > 0.05 & SDvND$p_val > 0.05)

mat <- matrix(c(a,b,c,d), byrow = T, nrow = 2, ncol = 2)

fisher.test(mat)
# Fisher's Exact Test for Count Data
# 
# data:  mat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  27.40339 38.10166
# sample estimates:
# odds ratio 
#   32.26221 

# Correlation of allele frequency differences between SD and ND vs SD and DIA 

SDvND <- read.table("./results/SDvND_freqDifs.txt", header=T)
SDvDIA <- read.table("./results/SDvCD_freqDifs.txt", header=T)

cor.test(SDvND$freqDif, SDvDIA$freqDif, method = "spearman")

# Correlations on chromosomes 1 - 3

#load in the maping information 
load("./data/Mapped_RAD_loci.Rdata")
chr1_3 <- chr1_ind_p | chr2_ind_p | chr3_ind_p

cor.test(SDvND$freqDif[chr1_3], SDvDIA$freqDif[chr1_3], method = "spearman")


############################
# Analysis of BSLMM output #
############################


setwd("/media/raglandlab/ExtraDrive4/Clines_3/results/output")

# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
model_est <- function(dat, filename){
  # pve -> proportion of phenotypic variance explained by the genotypes
  pve<-c("PVE", mean(dat$pve),quantile(dat$pve, probs=c(0.5,0.025,0.975)))
  
  #h -> hyperparameter for PVE
  h<-c("h", mean(dat$h), quantile(dat$h, probs=c(0.5,0.025,0.975)))
  
  # pge -> proportion of genetic variance explained by major effect loci
  pge<-c("PGE",mean(dat$pge),quantile(dat$pge, probs=c(0.5,0.025,0.975)))
  
  #rho -> hyperparameter for pge
  rho<-c("rho", mean(dat$rho), quantile(dat$rho, probs=c(0.5,0.025,0.975)))
  
  # pi -> proportion of variants with non-zero effects
  pi<-c("pi",mean(dat$pi),quantile(dat$pi, probs=c(0.5,0.025,0.975)))
  
  # n.gamma -> number of variants with major effect
  n.gamma<-c("n.gamma",mean(dat$n_gamma),quantile(dat$n_gamma, probs=c(0.5,0.025,0.975)))
  # ==============================================================================
  
  # get table of hyperparameters
  # ==============================================================================
  hyp.params.table<-as.data.frame(rbind(pve,h,pge,rho,pi,n.gamma),row.names=F)
  colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
  # show table
  hyp.params.table
  # write table to file
  write.table(hyp.params.table, file=filename, sep="\t", quote=F)
  # ==============================================================================
}

# plot traces and distributions of hyperparameters
#This will tell us whether our model was fit and converged
# ==============================================================================
#filename has to be pdf
hyper_trace <- function(dat, filename){
  pdf(file=filename, width=8.3,height=11.7)
  layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))
  # PVE - Percent variances in the phenotype explained by the genotype
  # ------------------------------------------------------------------------------
  plot(dat$pve, type="l", ylab="PVE", main="PVE - trace")
  hist(dat$pve, main="PVE - posterior distribution", xlab="PVE")
  plot(density(dat$pve), main="PVE - posterior distribution", xlab="PVE")
  # ------------------------------------------------------------------------------
  # PGE - proportion of PVE that is explained by loci of measurable effect
  # ------------------------------------------------------------------------------
  plot(dat$pge, type="l", ylab="PGE", main="PGE - trace")
  hist(dat$pge, main="PGE - posterior distribution", xlab="PGE")
  plot(density(dat$pge), main="PGE - posterior distribution", xlab="PGE")
  # ------------------------------------------------------------------------------
  # pi - proportion of Betas that are non zero (have an effect on the)
  # ------------------------------------------------------------------------------
  plot(dat$pi, type="l", ylab="pi", main="pi")
  hist(dat$pi, main="pi", xlab="pi")
  plot(density(dat$pi), main="pi", xlab="pi")
  # ------------------------------------------------------------------------------
  # h - Approximation to the proportion of phenotypic variance explained by genetic variantes (PVE)
  # ------------------------------------------------------------------------------
  plot(dat$h, type="l", ylab="h", main="h")
  hist(dat$h, main="h", xlab="h")
  plot(density(dat$h), main="h", xlab="h")
  # ------------------------------------------------------------------------------
  # rho - Approximation to the proportion of genetic variance explained by variants with a major effect (PGE)
  # ------------------------------------------------------------------------------
  plot(dat$rho, type="l", ylab="rho", main="rho")
  hist(dat$rho, main="rho", xlab="rho")
  plot(density(dat$rho), main="rho", xlab="rho")
  # ------------------------------------------------------------------------------
  # no. gamma - Number of loci of measureable effect (indicator variables that are 1)
  # ------------------------------------------------------------------------------
  plot(dat$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
  hist(dat$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
  plot(density(dat$pi), main="n_gamma - posterior distribution", xlab="n_gamma")
  # ------------------------------------------------------------------------------
  dev.off()
  
}

#Print the traces of the model runs and model estimates


#Low coverage individuals filtered, loci with missingness below 5% across all phenotpye groups - linear model
hyp.params<-read.table("SDCD_filtIndiv_filtSet_lm_1b5i.hyp.txt",header=T)
hyper_trace(dat = hyp.params, filename = "./param_traces_post/SDCD_filtIndiv_filtSet_lm_1b5i.hyp.pdf")
model_est(dat = hyp.params, filename = "./param_traces_post/SDCD_filtIndiv_filtSet_lm_1b5i.dsv")
#SD vs. ND 
hyp.params<-read.table("SDND_filtIndiv_filtSet_lm_1b5i.hyp.txt",header=T)
hyper_trace(dat = hyp.params, filename = "./param_traces_post/SDND_filtIndiv_filtSet_lm_1b5i.hyp.pdf")
model_est(dat = hyp.params, filename = "./param_traces_post/SDND_filtIndiv_filtSet_lm_1b5i.dsv")

#CD vs. ND
hyp.params<-read.table("CDND_filtIndiv_filtSet_lm_1b5i.hyp.txt",header=T)
hyper_trace(dat = hyp.params, filename = "./param_traces_post/CDND_filtIndiv_filtSet_lm_1b5i.hyp.pdf")
model_est(dat = hyp.params, filename = "./param_traces_post/CDND_filtIndiv_filtSet_lm_1b5i.dsv")



####violin plots of pve 
library(ggplot2)
CDvSD <- read.table("SDCD_filtIndiv_filtSet_lm_1b5i.hyp.txt", header=T)
NDvSD <- read.table("SDND_filtIndiv_filtSet_lm_1b5i.hyp.txt", header = T)
NDvCD <- read.table("CDND_filtIndiv_filtSet_lm_1b5i.hyp.txt", header = T)

comb <- rbind(NDvSD, CDvSD, NDvCD)

comb$trt <- c(rep("NDvSD", 125000), rep("CDvSD", 125000), rep("NDvCD", 125000))

mm_pve <- data.frame(mean = c(mean(comb$pve[comb$trt == "NDvSD"]), 
                              mean(comb$pve[comb$trt == "CDvSD"]),
                              mean(comb$pve[comb$trt == "NDvCD"])),
                     median = c(median(comb$pve[comb$trt == "NDvSD"]), 
                                median(comb$pve[comb$trt == "CDvSD"]),
                                median(comb$pve[comb$trt == "NDvCD"])),
                     trt = c("NDvSD", "CDvSD", "NDvCD"))

P <- ggplot(comb, aes(x=trt, y=pve, fill = trt)) + 
  geom_violin()
pdf("PVE_violin.pdf")
P + geom_boxplot(width=0.05, outlier.shape = NA, fill="black", coef = 0) + 
  geom_point(data = mm_pve, aes(x=trt, y=mean), color = "red")+
  geom_point(data = mm_pve, aes(x=trt, y=median), color = "white")+
  scale_x_discrete(limits=c("NDvSD", "CDvSD", "NDvCD")) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999","#999999","#999999"))+
  ylab("PVE") +
  xlab("Pairwise tests")+
  theme(legend.position="none")
dev.off()


mm_gamma <- data.frame(mean = c(mean(comb$n_gamma[comb$trt == "NDvSD"]), 
                              mean(comb$n_gamma[comb$trt == "CDvSD"]),
                              mean(comb$n_gamma[comb$trt == "NDvCD"])),
                     median = c(median(comb$n_gamma[comb$trt == "NDvSD"]), 
                                median(comb$n_gamma[comb$trt == "CDvSD"]),
                                median(comb$n_gamma[comb$trt == "NDvCD"])),
                     trt = c("NDvSD", "CDvSD", "NDvCD"))

### Violin plots of number of loci (n gamma)

P <- ggplot(comb, aes(x=trt, y=n_gamma, fill = trt)) + 
  geom_violin()
#pdf("PVE_violin.pdf")
P + geom_boxplot(width=0.05, outlier.shape = NA, fill="black", coef = 0) + 
  geom_point(data = mm_gamma, aes(x=trt, y=mean), color = "red")+
  geom_point(data = mm_gamma, aes(x=trt, y=median), color = "white")+
  scale_x_discrete(limits=c("NDvSD", "CDvSD", "NDvCD")) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999","#999999","#999999"))+
  ylab("PVE") +
  xlab("Pairwise tests")+
  theme(legend.position="none")
dev.off()


setwd("/media/raglandlab/ExtraDrive4/Clines_3/")

#Now for the polygenic score part 
#This script will be calcuated by firt determining which is the most common 
#allele in the SD group for each locus. E.g. if reverse phred scaled genotype 
#liklihood is > 1, then the alt is the more "SD" like allele, if < 1, then the "ref" is
#the more common allele. 

#The vcf file that is being loaded is already filtered by individual. 
library(vcfR)
Das <- read.vcfR("./data/snps.GATK.all.filtered.vcf", nrows = -1, skip = 0, cols = NULL, convertNA = TRUE, verbose = T)
PLs <- extract.gt(Das, element = "PL", mask = F, as.numeric = F, return.alleles = F, IDtoRowNames = TRUE, extract = T)
splitPls <- function(x){
  a <-as.numeric(unlist(strsplit(x,","))) #splits into vectro of three characters, then converts to numeric
  b <- 10^(-a/10)
  c <- 1/(sum(b))
  d <- c*b
  e <- sum(d*c(0,1,2))
  return(e) 
}
r <- nrow(PLs)
c <- ncol(PLs)
mat <- matrix(lapply(PLs, splitPls), nrow=r, ncol=c)
mat <- apply(mat,2,as.numeric)
colnames(mat)<-colnames(PLs)
rownames(mat)<-rownames(PLs)
SD_Ids <- read.table("./data/complete_shallowDiaIds.txt", header=F)
ND_Ids <- read.table("./data/complete_nonDiaIds.txt", header=F)
CD_Ids <- read.table("./data/complete_diaIds.txt", header=F)
bad_ind <- read.table("./data/remove_ind.txt", header=F) #individuals with very low read counts
badind <- colnames(mat)%in%bad_ind[,1]
mat <- mat[,!badind]
SDind <- colnames(mat)%in%SD_Ids[,1]
NDind <- colnames(mat)%in%ND_Ids[,1]
CDind <- colnames(mat)%in%CD_Ids[,1]
sd_mat <- mat[,SDind]
nd_mat <- mat[,NDind]
cd_mat <- mat[,CDind]
CDvSD <- read.table("./results/output/SDCD_filtIndiv_filtSet_lm_1b5i.param.txt", header=T)
NDvSD <- read.table("./results/output/SDND_filtIndiv_filtSet_lm_1b5i.param.txt", header = T)
NDvCD <- read.table("./results/output/CDND_filtIndiv_filtSet_lm_1b5i.param.txt", header = T)
ind <- NDvSD$gamma > 0.025
NDvSD_001_pos <- as.character(NDvSD$rs[ind])
ind <- CDvSD$gamma > 0.025
CDvSD_001_pos <- as.character(CDvSD$rs[ind])
ind <- NDvCD$gamma > 0.025
NDvCD_001_pos <- as.character(NDvCD$rs[ind])
pos_001 <- c(NDvSD_001_pos, CDvSD_001_pos, CDvSD_001_pos)

ind <- rownames(sd_mat) %in% pos_001
sd_mat_001 <-sd_mat[ind,]
ind <- rownames(nd_mat) %in% pos_001
nd_mat_001 <- nd_mat[ind,]
ind <- rownames(cd_mat) %in% pos_001
cd_mat_001 <- cd_mat[ind,]

freqEst <- function(f){
  if(sum(is.na(f)) == length(f)) { 
    nf = 1
    f = 0
  } else {
    f <- f[!is.na(f)]
    nf = length(f)
  }
  freq<-(sum(f)/(2*nf))
  return(freq)
}
indiv_score <- function(f){
  
}

sd_mat_001_freqs <- apply(sd_mat_001,1,freqEst)
nd_mat_001_freqs <- apply(nd_mat_001,1,freqEst)
cd_mat_001_freqs <- apply(cd_mat_001,1,freqEst)
ndcd_mat_001_freqs <- cbind(cd_mat_001_freqs, nd_mat_001_freqs)
ndcd_mat_001_freqs <- rowMeans(ndcd_mat_001_freqs)
comp_mat <- cbind(sd_mat_001_freqs, ndcd_mat_001_freqs)
comp_mat_2 <- comp_mat[,1] - comp_mat[,2]
pol_ind <- vector(length = length(comp_mat_2))
#true means that the alt allele is the most common, 
#false means that the reference allele is more common. 
for(i in 1:length(comp_mat_2)){
  if(comp_mat_2[i] > 0){
    pol_ind[i] <- TRUE
  }else{
    pol_ind[i] <- FALSE
  }
}
#change it up so the frequencies are those of the SD like allele
sd_poled <- matrix(nrow=nrow(sd_mat_001), ncol = ncol(sd_mat_001))
for(i in 1:length(pol_ind)){
  for(j in 1:ncol(sd_poled))
    if(pol_ind[i]==TRUE){
      sd_poled[i,j] <- sd_mat_001[i,j]
    }else if(pol_ind[i]==FALSE){
      sd_poled[i,j] <- (2-sd_mat_001[i,j])
    }
}
sd_scored<- apply(sd_poled,2,freqEst)

cd_poled <- matrix(nrow=nrow(cd_mat_001), ncol = ncol(cd_mat_001))
for(i in 1:length(pol_ind)){
  for(j in 1:ncol(cd_poled))
    if(pol_ind[i]==TRUE){
      cd_poled[i,j] <- cd_mat_001[i,j]
    }else if(pol_ind[i]==FALSE){
      cd_poled[i,j] <- (2-cd_mat_001[i,j])
    }
}
cd_scored <- apply(cd_poled,2,freqEst)

nd_poled <- matrix(nrow=nrow(nd_mat_001), ncol = ncol(nd_mat_001))
for(i in 1:length(pol_ind)){
  for(j in 1:ncol(nd_poled))
    if(pol_ind[i]==TRUE){
      nd_poled[i,j] <- nd_mat_001[i,j]
    }else if(pol_ind[i]==FALSE){
      nd_poled[i,j] <- (2-nd_mat_001[i,j])
    }
}
nd_scored <- apply(nd_poled,2,freqEst)

pdf("polygenescore_bslmm001_gamma.pdf")
plot(density(cd_scored, bw=0.015), xlim = c(0.3,0.6), ylim = c(0,17), main=NA, col = "#1b9e77", lty=4,lwd=2)
lines(density(nd_scored, bw=0.015), col = "#d95f02", lty=1,lwd=2)
lines(density(sd_scored, bw=0.015), col="#7570b3", lty=2, lwd=2)
rug(cd_scored, col="#1b9e77", lwd = 1)
rug(nd_scored, col="#d95f02", lwd = 1)
rug(sd_scored, col = "#7570b3", lwd = 1)
legend("topright", legend=c("non diapause", "diapause", "shallow"), col=c("#d95f02", "#1b9e77", "#7570b3"), lty=c(1,4,2), lwd=2)
dev.off()

#####K-S test among polygenic scores to test for differences in the distribution 

#SD vs. Dia 
ks.test(sd_scored, cd_scored)
# Two-sample Kolmogorov-Smirnov test
# 
# data:  sd_scored and cd_scored
# D = 0.83448, p-value < 2.2e-16
# alternative hypothesis: two-sided

#SD vs. ND
ks.test(sd_scored, nd_scored)
# Two-sample Kolmogorov-Smirnov test
# 
# data:  sd_scored and nd_scored
# D = 0.92034, p-value < 2.2e-16
# alternative hypothesis: two-sided

#Dia vs. ND
ks.test(cd_scored, nd_scored)
# Two-sample Kolmogorov-Smirnov test
# 
# data:  cd_scored and nd_scored
# D = 0.10336, p-value = 0.8632
# alternative hypothesis: two-sided

######## Polygenic score for all ~7k loci in the shared data set. 
#sd_mat, nd_mat, and cd_mat are the dataframes with the mean genotype freqs for all loci

#Filter down to ~7k loci 
shared_loci <- read.table("./data/position_fix.txt")
sd_mat <- subset(sd_mat, rownames(sd_mat) %in% shared_loci$V1)
nd_mat <- subset(nd_mat, rownames(nd_mat) %in% shared_loci$V1)
cd_mat <- subset(cd_mat, rownames(cd_mat) %in% shared_loci$V1)


sd_mat_freqs <- apply(sd_mat,1,freqEst)
nd_mat_freqs <- apply(nd_mat,1,freqEst)
cd_mat_freqs <- apply(cd_mat,1,freqEst)
ndcd_mat_freqs <- cbind(cd_mat_freqs, nd_mat_freqs)
ndcd_mat_freqs <- rowMeans(ndcd_mat_freqs)
comp_mat <- cbind(sd_mat_freqs, ndcd_mat_freqs)
comp_mat_2 <- comp_mat[,1] - comp_mat[,2]
pol_ind <- vector(length = length(comp_mat_2))
#true means that the alt allele is the most common, 
#false means that the reference allele is more common. 
for(i in 1:length(comp_mat_2)){
  if(comp_mat_2[i] > 0){
    pol_ind[i] <- TRUE
  }else{
    pol_ind[i] <- FALSE
  }
}
#change it up so the frequencies are those of the SD like allele
sd_poled <- matrix(nrow=nrow(sd_mat), ncol = ncol(sd_mat))
for(i in 1:length(pol_ind)){
  for(j in 1:ncol(sd_poled))
    if(pol_ind[i]==TRUE){
      sd_poled[i,j] <- sd_mat[i,j]
    }else if(pol_ind[i]==FALSE){
      sd_poled[i,j] <- (2-sd_mat[i,j])
    }
}
sd_scored<- apply(sd_poled,2,freqEst)

cd_poled <- matrix(nrow=nrow(cd_mat), ncol = ncol(cd_mat))
for(i in 1:length(pol_ind)){
  for(j in 1:ncol(cd_poled))
    if(pol_ind[i]==TRUE){
      cd_poled[i,j] <- cd_mat[i,j]
    }else if(pol_ind[i]==FALSE){
      cd_poled[i,j] <- (2-cd_mat[i,j])
    }
}
cd_scored <- apply(cd_poled,2,freqEst)

nd_poled <- matrix(nrow=nrow(nd_mat), ncol = ncol(nd_mat))
for(i in 1:length(pol_ind)){
  for(j in 1:ncol(nd_poled))
    if(pol_ind[i]==TRUE){
      nd_poled[i,j] <- nd_mat[i,j]
    }else if(pol_ind[i]==FALSE){
      nd_poled[i,j] <- (2-nd_mat[i,j])
    }
}
nd_scored <- apply(nd_poled,2,freqEst)


#We see the same pattern for full data set, so we will use this figure because it is more informative. 
pdf("./results/DiaClass_polygenicScoreAllLoci.pdf")
plot(density(cd_scored, bw=0.015), xlim = c(0.3,0.6), ylim = c(0,17), main=NA, col = "#1b9e77", lty=4,lwd=2)
lines(density(nd_scored, bw=0.015), col = "#d95f02", lty=1,lwd=2)
lines(density(sd_scored, bw=0.015), col="#7570b3", lty=2, lwd=2)
rug(cd_scored, col="#1b9e77", lwd = 1)
rug(nd_scored, col="#d95f02", lwd = 1)
rug(sd_scored, col = "#7570b3", lwd = 1)
legend("topright", legend=c("non diapause", "diapause", "shallow"), col=c("#d95f02", "#1b9e77", "#7570b3"), lty=c(1,4,2), lwd=2)
dev.off()


