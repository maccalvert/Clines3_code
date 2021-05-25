# This set of scripts is to calculate the frequencies of allele associated with diapause intensity and alleles associated with 
# diapause termination at the four geographic locations in this study. 

library(doParallel)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(scales)
library(VennDiagram)

setwd("/media/raglandlab/ExtraDrive4/Clines_3")
source("/media/raglandlab/ExtraDrive4/Clines_3/src/clines3_functions.R")

Gra_haw <- read.table("./data/haw7_genosAll_polarized.txt")
Gra_apple <- read.table("./data/apple7_genosAll_polarized.txt")
Fen_haw <- read.table("./data/FenHaw_genosAllChr_polarized.txt", header=T)
Fen_apple <- read.table("./data/FenApple_genosAllChr_polarized.txt", header=T)
Dow_haw <- read.table("./data/DowHaw_genosAllChr_polarized.txt", header=T)
Dow_apple <- read.table("./data/DowApple_genosAllChr_polarized.txt", header=T)
Urb_haw <- read.table("./data/UrbHaw_genosAllChr_polarized.txt", header=T)
Urb_apple <- read.table("./data/UrbApple_genosAllChr_polarized.txt", header=T)

cl <- makeCluster(10)
registerDoParallel(cl)
nloci=7265
freqEst_GraHaw<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Gra_haw[i,]))
}
freqEst_GraApple<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Gra_apple[i,]))
}
freqEst_FenHaw<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Fen_haw[i,]))
}
freqEst_FenApple<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Fen_apple[i,]))
}
freqEst_DowHaw<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Dow_haw[i,]))
}
freqEst_DowApple<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Dow_apple[i,]))
}
freqEst_UrbHaw<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Urb_haw[i,]))
}
freqEst_UrbApple<-foreach(i=1:nloci, .combine=c) %dopar% {
  freqEst(unlist(Urb_apple[i,]))
}
stopCluster(cl)
# save.image("FreqEst_clines.Rdata")
# load("FreqEst_clines.Rdata")

load("./data/Mapped_RAD_loci.Rdata")

AvgEcl <- read.table("./results/avg_early_late_freqDifs.txt", header=T)
#CDNDvSD <- read.table("SDvCDND_freqDifs.txt", header=T)
DiaInt <- read.table("./results/SDvCDND_freqDifs.txt", header=T)
colnames(DiaInt) <- c("freqDif", "pval")
rownames(DiaInt) <- rownames(AvgEcl)
Ecl_loci <- AvgEcl[AvgEcl$p_val < 0.05,]
Dia_loci <- DiaInt[DiaInt$pval < 0.05,]

#This section polarizes loci so that we are counting the allele most associated with the shallow diapause phenotype, or in the 
# the case of the diapause termination study, we are counting the allele most associated with late diapause phenotype. 

SD_pol <- vector(length=nrow(DiaInt))
for(i in 1:nrow(DiaInt)){
  if(DiaInt$freqDif[i] < 0){
    SD_pol[i] <- T
  }else if(DiaInt$freqDif[i] > 0){
    SD_pol[i] <- F
  }
}

late_pol <- vector(length=nrow(AvgEcl))
for(i in 1:nrow(AvgEcl)){
  if(AvgEcl$freqDif[i] > 0){
    late_pol[i] <- T
  }else if(AvgEcl$freqDif[i] < 0){
    late_pol[i] <- F
  }
}

#SD polarized vectors 
freqEst_GraHaw_pol <- pol_func(vec=freqEst_GraHaw, pol=unlist(SD_pol))
freqEst_GraApple_pol <- pol_func(vec=freqEst_GraApple, pol=unlist(SD_pol))
freqEst_FenHaw_pol <- pol_func(vec=freqEst_FenHaw, pol=unlist(SD_pol))
freqEst_FenApple_pol <- pol_func(vec=freqEst_FenApple, pol=unlist(SD_pol))
freqEst_DowHaw_pol <- pol_func(vec=freqEst_DowHaw, pol=unlist(SD_pol))
freqEst_DowApple_pol <- pol_func(vec=freqEst_DowApple, pol=unlist(SD_pol))
freqEst_UrbHaw_pol <- pol_func(vec=freqEst_UrbHaw, pol=unlist(SD_pol))
freqEst_UrbApple_pol <- pol_func(vec=freqEst_UrbApple, pol=unlist(SD_pol))

#late polarized vectors 
freqEst_GraHaw_late_pol <- pol_func(vec=freqEst_GraHaw, pol=unlist(late_pol))
freqEst_GraApple_late_pol <- pol_func(vec=freqEst_GraApple, pol=unlist(late_pol))
freqEst_FenHaw_late_pol <- pol_func(vec=freqEst_FenHaw, pol=unlist(late_pol))
freqEst_FenApple_late_pol <- pol_func(vec=freqEst_FenApple, pol=unlist(late_pol))
freqEst_DowHaw_late_pol <- pol_func(vec=freqEst_DowHaw, pol=unlist(late_pol))
freqEst_DowApple_late_pol <- pol_func(vec=freqEst_DowApple, pol=unlist(late_pol))
freqEst_UrbHaw_late_pol <- pol_func(vec=freqEst_UrbHaw, pol=unlist(late_pol))
freqEst_UrbApple_late_pol <- pol_func(vec=freqEst_UrbApple, pol=unlist(late_pol))

#Venn diagrams of overlapping loci.

#Chromosome 1 

ecl_ind <- chr1_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr1_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

mypallete <- viridis_pal(option = "D")(1000)

snps_1 <- rownames(DiaInt)[dia_ind]

snps_2 <- rownames(AvgEcl)[ecl_ind]

pdf("./results/ven_diaecl_Chr1.pdf")
grid.draw(venn.diagram(x=list(snps_1, snps_2), 
                       category.names=c("Diapasue Intesntiy" , "Eclosion timing"),
                       fill = c(mypallete[100], mypallete[400]),
                       filename=NULL))
dev.off()

#Chromosome 2 

ecl_ind <- chr2_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr2_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

mypallete <- viridis_pal(option = "D")(1000)

snps_1 <- rownames(DiaInt)[dia_ind]
snps_2 <- rownames(AvgEcl)[ecl_ind]

pdf("./results/ven_diaecl_Chr2.pdf")
grid.draw(venn.diagram(x=list(snps_1, snps_2), 
                       category.names=c("Diapasue Intesntiy" , "Eclosion timing"),
                       fill = c(mypallete[100], mypallete[400]),
                       filename=NULL))
dev.off()

#Chromosome 3

ecl_ind <- chr3_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr3_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

mypallete <- viridis_pal(option = "D")(1000)

snps_1 <- rownames(DiaInt)[dia_ind]
snps_2 <- rownames(AvgEcl)[ecl_ind]

pdf("./results/ven_diaecl_Chr3.pdf")
grid.draw(venn.diagram(x=list(snps_1, snps_2), 
                       category.names=c("Diapasue Intesntiy" , "Eclosion timing"),
                       fill = c(mypallete[100], mypallete[400]),
                       filename=NULL))
dev.off()

#Chromosome 4

ecl_ind <- chr4_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr4_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

mypallete <- viridis_pal(option = "D")(1000)

snps_1 <- rownames(DiaInt)[dia_ind]
snps_2 <- rownames(AvgEcl)[ecl_ind]

pdf("./results/ven_diaecl_Chr4.pdf")
grid.draw(venn.diagram(x=list(snps_1, snps_2), 
                       category.names=c("Diapasue Intesntiy" , "Eclosion timing"),
                       fill = c(mypallete[100], mypallete[400]),
                       filename=NULL))
dev.off()


#Chromosome 5

ecl_ind <- chr5_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr5_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

mypallete <- viridis_pal(option = "D")(1000)

snps_1 <- rownames(DiaInt)[dia_ind]
snps_2 <- rownames(AvgEcl)[ecl_ind]

pdf("./results/ven_diaecl_Chr5.pdf")
grid.draw(venn.diagram(x=list(snps_1, snps_2), 
                       category.names=c("Diapasue Intesntiy" , "Eclosion timing"),
                       fill = c(mypallete[100], mypallete[400]),
                       filename=NULL))
dev.off()


#Polarized allele frequency clines. 

#Chromosome 1

#SD polarized allele freq figure - plieotropic and SD only 

ecl_ind <- chr1_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr1_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

Dia_means <- c(mean(freqEst_GraHaw_pol[dia_ind]), mean(freqEst_GraApple_pol[dia_ind]),
                    mean(freqEst_FenHaw_pol[dia_ind]), mean(freqEst_FenApple_pol[dia_ind]),
                    mean(freqEst_DowHaw_pol[dia_ind]), mean(freqEst_DowApple_pol[dia_ind]),
                    mean(freqEst_UrbHaw_pol[dia_ind]), mean(freqEst_UrbApple_pol[dia_ind]))
Dia_CI <- rbind(CI_95(freqEst_GraHaw_pol[dia_ind], n=length(freqEst_GraHaw_pol[dia_ind])), CI_95(freqEst_GraApple_pol[dia_ind], n=length(freqEst_GraApple_pol[dia_ind])),
                     CI_95(freqEst_FenHaw_pol[dia_ind], n=length(freqEst_FenHaw_pol[dia_ind])), CI_95(freqEst_FenApple_pol[dia_ind], n=length(freqEst_FenApple_pol[dia_ind])),
                     CI_95(freqEst_DowHaw_pol[dia_ind], n=length(freqEst_DowHaw_pol[dia_ind])), CI_95(freqEst_DowApple_pol[dia_ind], n=length(freqEst_DowApple_pol[dia_ind])),
                     CI_95(freqEst_UrbHaw_pol[dia_ind], n=length(freqEst_UrbHaw_pol[dia_ind])), CI_95(freqEst_UrbApple_pol[dia_ind], n=length(freqEst_UrbApple_pol[dia_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Dia_means, Dia_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_dia_Chr1.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()

#late polarized allele freq figure. - plieotropic and late only 

Ecl_means <- c(mean(freqEst_GraHaw_late_pol[ecl_ind]), mean(freqEst_GraApple_late_pol[ecl_ind]),
                    mean(freqEst_FenHaw_late_pol[ecl_ind]), mean(freqEst_FenApple_late_pol[ecl_ind]),
                    mean(freqEst_DowHaw_late_pol[ecl_ind]), mean(freqEst_DowApple_late_pol[ecl_ind]),
                    mean(freqEst_UrbHaw_late_pol[ecl_ind]), mean(freqEst_UrbApple_late_pol[ecl_ind]))
Ecl_CI <- rbind(CI_95(freqEst_GraHaw_late_pol[ecl_ind], n=length(freqEst_GraHaw_late_pol[ecl_ind])), CI_95(freqEst_GraApple_late_pol[ecl_ind], n=length(freqEst_GraApple_late_pol[ecl_ind])),
                     CI_95(freqEst_FenHaw_late_pol[ecl_ind], n=length(freqEst_FenHaw_late_pol[ecl_ind])), CI_95(freqEst_FenApple_late_pol[ecl_ind], n=length(freqEst_FenApple_late_pol[ecl_ind])),
                     CI_95(freqEst_DowHaw_late_pol[ecl_ind], n=length(freqEst_DowHaw_late_pol[ecl_ind])), CI_95(freqEst_DowApple_late_pol[ecl_ind], n=length(freqEst_DowApple_late_pol[ecl_ind])),
                     CI_95(freqEst_UrbHaw_late_pol[ecl_ind], n=length(freqEst_UrbHaw_late_pol[ecl_ind])), CI_95(freqEst_UrbApple_late_pol[ecl_ind], n=length(freqEst_UrbApple_late_pol[ecl_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Ecl_means, Ecl_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_ecl_Chr1.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()


#Chromosome 2

#SD polarized allele freq figure - plieotropic and SD only 

ecl_ind <- chr2_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr2_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

Dia_means <- c(mean(freqEst_GraHaw_pol[dia_ind]), mean(freqEst_GraApple_pol[dia_ind]),
               mean(freqEst_FenHaw_pol[dia_ind]), mean(freqEst_FenApple_pol[dia_ind]),
               mean(freqEst_DowHaw_pol[dia_ind]), mean(freqEst_DowApple_pol[dia_ind]),
               mean(freqEst_UrbHaw_pol[dia_ind]), mean(freqEst_UrbApple_pol[dia_ind]))
Dia_CI <- rbind(CI_95(freqEst_GraHaw_pol[dia_ind], n=length(freqEst_GraHaw_pol[dia_ind])), CI_95(freqEst_GraApple_pol[dia_ind], n=length(freqEst_GraApple_pol[dia_ind])),
                CI_95(freqEst_FenHaw_pol[dia_ind], n=length(freqEst_FenHaw_pol[dia_ind])), CI_95(freqEst_FenApple_pol[dia_ind], n=length(freqEst_FenApple_pol[dia_ind])),
                CI_95(freqEst_DowHaw_pol[dia_ind], n=length(freqEst_DowHaw_pol[dia_ind])), CI_95(freqEst_DowApple_pol[dia_ind], n=length(freqEst_DowApple_pol[dia_ind])),
                CI_95(freqEst_UrbHaw_pol[dia_ind], n=length(freqEst_UrbHaw_pol[dia_ind])), CI_95(freqEst_UrbApple_pol[dia_ind], n=length(freqEst_UrbApple_pol[dia_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Dia_means, Dia_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_dia_Chr2.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()

#late polarized allele freq figure. - plieotropic and late only 

Ecl_means <- c(mean(freqEst_GraHaw_late_pol[ecl_ind]), mean(freqEst_GraApple_late_pol[ecl_ind]),
               mean(freqEst_FenHaw_late_pol[ecl_ind]), mean(freqEst_FenApple_late_pol[ecl_ind]),
               mean(freqEst_DowHaw_late_pol[ecl_ind]), mean(freqEst_DowApple_late_pol[ecl_ind]),
               mean(freqEst_UrbHaw_late_pol[ecl_ind]), mean(freqEst_UrbApple_late_pol[ecl_ind]))
Ecl_CI <- rbind(CI_95(freqEst_GraHaw_late_pol[ecl_ind], n=length(freqEst_GraHaw_late_pol[ecl_ind])), CI_95(freqEst_GraApple_late_pol[ecl_ind], n=length(freqEst_GraApple_late_pol[ecl_ind])),
                CI_95(freqEst_FenHaw_late_pol[ecl_ind], n=length(freqEst_FenHaw_late_pol[ecl_ind])), CI_95(freqEst_FenApple_late_pol[ecl_ind], n=length(freqEst_FenApple_late_pol[ecl_ind])),
                CI_95(freqEst_DowHaw_late_pol[ecl_ind], n=length(freqEst_DowHaw_late_pol[ecl_ind])), CI_95(freqEst_DowApple_late_pol[ecl_ind], n=length(freqEst_DowApple_late_pol[ecl_ind])),
                CI_95(freqEst_UrbHaw_late_pol[ecl_ind], n=length(freqEst_UrbHaw_late_pol[ecl_ind])), CI_95(freqEst_UrbApple_late_pol[ecl_ind], n=length(freqEst_UrbApple_late_pol[ecl_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Ecl_means, Ecl_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_ecl_Chr2.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()


#Chromosome 3

#SD polarized allele freq figure - plieotropic and SD only 

ecl_ind <- chr3_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr3_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

Dia_means <- c(mean(freqEst_GraHaw_pol[dia_ind]), mean(freqEst_GraApple_pol[dia_ind]),
               mean(freqEst_FenHaw_pol[dia_ind]), mean(freqEst_FenApple_pol[dia_ind]),
               mean(freqEst_DowHaw_pol[dia_ind]), mean(freqEst_DowApple_pol[dia_ind]),
               mean(freqEst_UrbHaw_pol[dia_ind]), mean(freqEst_UrbApple_pol[dia_ind]))
Dia_CI <- rbind(CI_95(freqEst_GraHaw_pol[dia_ind], n=length(freqEst_GraHaw_pol[dia_ind])), CI_95(freqEst_GraApple_pol[dia_ind], n=length(freqEst_GraApple_pol[dia_ind])),
                CI_95(freqEst_FenHaw_pol[dia_ind], n=length(freqEst_FenHaw_pol[dia_ind])), CI_95(freqEst_FenApple_pol[dia_ind], n=length(freqEst_FenApple_pol[dia_ind])),
                CI_95(freqEst_DowHaw_pol[dia_ind], n=length(freqEst_DowHaw_pol[dia_ind])), CI_95(freqEst_DowApple_pol[dia_ind], n=length(freqEst_DowApple_pol[dia_ind])),
                CI_95(freqEst_UrbHaw_pol[dia_ind], n=length(freqEst_UrbHaw_pol[dia_ind])), CI_95(freqEst_UrbApple_pol[dia_ind], n=length(freqEst_UrbApple_pol[dia_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Dia_means, Dia_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_dia_Chr3.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()

#late polarized allele freq figure. - plieotropic and late only 

Ecl_means <- c(mean(freqEst_GraHaw_late_pol[ecl_ind]), mean(freqEst_GraApple_late_pol[ecl_ind]),
               mean(freqEst_FenHaw_late_pol[ecl_ind]), mean(freqEst_FenApple_late_pol[ecl_ind]),
               mean(freqEst_DowHaw_late_pol[ecl_ind]), mean(freqEst_DowApple_late_pol[ecl_ind]),
               mean(freqEst_UrbHaw_late_pol[ecl_ind]), mean(freqEst_UrbApple_late_pol[ecl_ind]))
Ecl_CI <- rbind(CI_95(freqEst_GraHaw_late_pol[ecl_ind], n=length(freqEst_GraHaw_late_pol[ecl_ind])), CI_95(freqEst_GraApple_late_pol[ecl_ind], n=length(freqEst_GraApple_late_pol[ecl_ind])),
                CI_95(freqEst_FenHaw_late_pol[ecl_ind], n=length(freqEst_FenHaw_late_pol[ecl_ind])), CI_95(freqEst_FenApple_late_pol[ecl_ind], n=length(freqEst_FenApple_late_pol[ecl_ind])),
                CI_95(freqEst_DowHaw_late_pol[ecl_ind], n=length(freqEst_DowHaw_late_pol[ecl_ind])), CI_95(freqEst_DowApple_late_pol[ecl_ind], n=length(freqEst_DowApple_late_pol[ecl_ind])),
                CI_95(freqEst_UrbHaw_late_pol[ecl_ind], n=length(freqEst_UrbHaw_late_pol[ecl_ind])), CI_95(freqEst_UrbApple_late_pol[ecl_ind], n=length(freqEst_UrbApple_late_pol[ecl_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Ecl_means, Ecl_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_ecl_Chr3.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()


#Chromosome 4

#SD polarized allele freq figure - plieotropic and SD only 

ecl_ind <- chr4_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr4_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

Dia_means <- c(mean(freqEst_GraHaw_pol[dia_ind]), mean(freqEst_GraApple_pol[dia_ind]),
               mean(freqEst_FenHaw_pol[dia_ind]), mean(freqEst_FenApple_pol[dia_ind]),
               mean(freqEst_DowHaw_pol[dia_ind]), mean(freqEst_DowApple_pol[dia_ind]),
               mean(freqEst_UrbHaw_pol[dia_ind]), mean(freqEst_UrbApple_pol[dia_ind]))
Dia_CI <- rbind(CI_95(freqEst_GraHaw_pol[dia_ind], n=length(freqEst_GraHaw_pol[dia_ind])), CI_95(freqEst_GraApple_pol[dia_ind], n=length(freqEst_GraApple_pol[dia_ind])),
                CI_95(freqEst_FenHaw_pol[dia_ind], n=length(freqEst_FenHaw_pol[dia_ind])), CI_95(freqEst_FenApple_pol[dia_ind], n=length(freqEst_FenApple_pol[dia_ind])),
                CI_95(freqEst_DowHaw_pol[dia_ind], n=length(freqEst_DowHaw_pol[dia_ind])), CI_95(freqEst_DowApple_pol[dia_ind], n=length(freqEst_DowApple_pol[dia_ind])),
                CI_95(freqEst_UrbHaw_pol[dia_ind], n=length(freqEst_UrbHaw_pol[dia_ind])), CI_95(freqEst_UrbApple_pol[dia_ind], n=length(freqEst_UrbApple_pol[dia_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Dia_means, Dia_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_dia_Chr4.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()

#late polarized allele freq figure. - plieotropic and late only 

Ecl_means <- c(mean(freqEst_GraHaw_late_pol[ecl_ind]), mean(freqEst_GraApple_late_pol[ecl_ind]),
               mean(freqEst_FenHaw_late_pol[ecl_ind]), mean(freqEst_FenApple_late_pol[ecl_ind]),
               mean(freqEst_DowHaw_late_pol[ecl_ind]), mean(freqEst_DowApple_late_pol[ecl_ind]),
               mean(freqEst_UrbHaw_late_pol[ecl_ind]), mean(freqEst_UrbApple_late_pol[ecl_ind]))
Ecl_CI <- rbind(CI_95(freqEst_GraHaw_late_pol[ecl_ind], n=length(freqEst_GraHaw_late_pol[ecl_ind])), CI_95(freqEst_GraApple_late_pol[ecl_ind], n=length(freqEst_GraApple_late_pol[ecl_ind])),
                CI_95(freqEst_FenHaw_late_pol[ecl_ind], n=length(freqEst_FenHaw_late_pol[ecl_ind])), CI_95(freqEst_FenApple_late_pol[ecl_ind], n=length(freqEst_FenApple_late_pol[ecl_ind])),
                CI_95(freqEst_DowHaw_late_pol[ecl_ind], n=length(freqEst_DowHaw_late_pol[ecl_ind])), CI_95(freqEst_DowApple_late_pol[ecl_ind], n=length(freqEst_DowApple_late_pol[ecl_ind])),
                CI_95(freqEst_UrbHaw_late_pol[ecl_ind], n=length(freqEst_UrbHaw_late_pol[ecl_ind])), CI_95(freqEst_UrbApple_late_pol[ecl_ind], n=length(freqEst_UrbApple_late_pol[ecl_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Ecl_means, Ecl_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_ecl_Chr4.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()


#Chromosome 5

#SD polarized allele freq figure - plieotropic and SD only 

ecl_ind <- chr5_ind_p & rownames(Gra_haw) %in% rownames(Ecl_loci)
dia_ind <- chr5_ind_p & rownames(Gra_haw) %in% rownames(Dia_loci)

Dia_means <- c(mean(freqEst_GraHaw_pol[dia_ind]), mean(freqEst_GraApple_pol[dia_ind]),
               mean(freqEst_FenHaw_pol[dia_ind]), mean(freqEst_FenApple_pol[dia_ind]),
               mean(freqEst_DowHaw_pol[dia_ind]), mean(freqEst_DowApple_pol[dia_ind]),
               mean(freqEst_UrbHaw_pol[dia_ind]), mean(freqEst_UrbApple_pol[dia_ind]))
Dia_CI <- rbind(CI_95(freqEst_GraHaw_pol[dia_ind], n=length(freqEst_GraHaw_pol[dia_ind])), CI_95(freqEst_GraApple_pol[dia_ind], n=length(freqEst_GraApple_pol[dia_ind])),
                CI_95(freqEst_FenHaw_pol[dia_ind], n=length(freqEst_FenHaw_pol[dia_ind])), CI_95(freqEst_FenApple_pol[dia_ind], n=length(freqEst_FenApple_pol[dia_ind])),
                CI_95(freqEst_DowHaw_pol[dia_ind], n=length(freqEst_DowHaw_pol[dia_ind])), CI_95(freqEst_DowApple_pol[dia_ind], n=length(freqEst_DowApple_pol[dia_ind])),
                CI_95(freqEst_UrbHaw_pol[dia_ind], n=length(freqEst_UrbHaw_pol[dia_ind])), CI_95(freqEst_UrbApple_pol[dia_ind], n=length(freqEst_UrbApple_pol[dia_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Dia_means, Dia_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_dia_Chr5.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()

#late polarized allele freq figure. - plieotropic and late only 

Ecl_means <- c(mean(freqEst_GraHaw_late_pol[ecl_ind]), mean(freqEst_GraApple_late_pol[ecl_ind]),
               mean(freqEst_FenHaw_late_pol[ecl_ind]), mean(freqEst_FenApple_late_pol[ecl_ind]),
               mean(freqEst_DowHaw_late_pol[ecl_ind]), mean(freqEst_DowApple_late_pol[ecl_ind]),
               mean(freqEst_UrbHaw_late_pol[ecl_ind]), mean(freqEst_UrbApple_late_pol[ecl_ind]))
Ecl_CI <- rbind(CI_95(freqEst_GraHaw_late_pol[ecl_ind], n=length(freqEst_GraHaw_late_pol[ecl_ind])), CI_95(freqEst_GraApple_late_pol[ecl_ind], n=length(freqEst_GraApple_late_pol[ecl_ind])),
                CI_95(freqEst_FenHaw_late_pol[ecl_ind], n=length(freqEst_FenHaw_late_pol[ecl_ind])), CI_95(freqEst_FenApple_late_pol[ecl_ind], n=length(freqEst_FenApple_late_pol[ecl_ind])),
                CI_95(freqEst_DowHaw_late_pol[ecl_ind], n=length(freqEst_DowHaw_late_pol[ecl_ind])), CI_95(freqEst_DowApple_late_pol[ecl_ind], n=length(freqEst_DowApple_late_pol[ecl_ind])),
                CI_95(freqEst_UrbHaw_late_pol[ecl_ind], n=length(freqEst_UrbHaw_late_pol[ecl_ind])), CI_95(freqEst_UrbApple_late_pol[ecl_ind], n=length(freqEst_UrbApple_late_pol[ecl_ind])))
site <- c(rep(c("Grant"),2), rep(c("Fennville"),2), rep(c("Dowagiac"),2), rep("Urbana",2))
race <- c(rep(c("Haw", "Apple"),8))
#phenotype <- c(rep("Eclosion", 8), rep("Diapause",8))
mypalette <- brewer.pal(name = "PiYG", n=10)

freq_df <- as.data.frame(cbind(Ecl_means, Ecl_CI, site,race))
colnames(freq_df) <- c("AlleleFreq", "CI2.5", "CI97.5", "Sites", "Race")
freq_df$AlleleFreq <- as.numeric(as.character(freq_df$AlleleFreq))
freq_df$CI2.5 <- as.numeric(as.character(freq_df$CI2.5))
freq_df$CI97.5 <- as.numeric(as.character(freq_df$CI97.5))
pdf("./results/meanAlleleFreq_ecl_Chr5.pdf")
p<- ggplot(data=freq_df, aes(x=Sites, y=AlleleFreq, color=Race, group=Race)) +
  scale_color_manual(values=c(mypalette[2], mypalette[8]))+
  geom_line(position = position_dodge(width=0.2), size=1.0) +
  geom_point(position = position_dodge(width=0.2), size=2.0) + 
  xlim("Urbana", "Dowagiac", "Fennville","Grant")+
  ylim(0.2, 0.75) +
  theme_classic() +
  #theme(axis.line = element_line(size=2.0), axis.text = element_text(size=20.0))+
  geom_errorbar(aes(ymin=CI2.5, ymax=CI97.5), width=.2, position=position_dodge(width=0.2), size=1.0)+
  coord_fixed(ratio = 6)
p
dev.off()
