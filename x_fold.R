### x_fold test to determine if there is more sign overlap between alleles
# frequencies in diapause and eclosion than would be expected by chance
#v1 
#MBC on 9/1/2018
#x-fold functions largely derived from Chaturvedi et al. 2018 

#Revised for online repository 
#MBC on 6/1/2020

#First need to load in genotypes and calculate permuated whole genotypes. 

library(iterators)
library(doParallel)
library(foreach)

setwd("/media/raglandlab/ExtraDrive4/Clines_3")
source("/media/raglandlab/ExtraDrive4/Clines_3/src/clines3_functions.R")

#load in the genotypes that have been polarized to Grant - Haw population, the same poiarization
#scheme as the rest of the paper. 

CDND <- read.table("./data/CDND_genosAll_polarized.txt", header=T)
SD <- read.table("./data/SD_genosAll_polarized.txt", header=T)
appleEarly <- read.table("./data/appleEarly_genosAll_polarized.txt", header=T)
appleLate <- read.table("./data/appleLate_genosAll_polarized.txt", header=T)
HawEarly <- read.table("./data/hawEarly_genosAll_polarized.txt", header=T)
HawLate <- read.table("./data/hawLate_genosAll_polarized.txt", header=T)

# Next, to the x-fold function is written to accept the observed allele frequency differences as table with rows representing 
# loci and columns representing allele frequency differences for each experiment. There also needs to be a column for chromosome 
# and LD group. 
#These respresent the "observed values" which will be testing against the null produced from 
#permuted allele frequency differences. 
Dia_freqDifs <- read.table("./results/SDvCDND_freqDifs.txt", header=T)
Ecl_freqDifs <- read.table("./results/avg_early_late_freqDifs.txt", header=T)

#bind them up and make empty columns for chromosome and LD group
load("./data/Mapped_RAD_loci.Rdata")

Alleles <- data.frame(chr = rep(NA, nrow(Ecl_freqDifs)), LDgr = rep(NA, nrow(Ecl_freqDifs)), 
                      "Dia" = Dia_freqDifs$freqDif, "Ecl" = Ecl_freqDifs$freqDif)
Alleles$chr[chr1_ind_p] <- 1
Alleles$chr[chr2_ind_p] <- 2
Alleles$chr[chr3_ind_p] <- 3
Alleles$chr[chr4_ind_p] <- 4
Alleles$chr[chr5_ind_p] <- 5
Alleles$LDgr[all_LDL_p] <- "L"
Alleles$LDgr[all_LDM_p] <- "M"
Alleles$LDgr[all_LDH_p] <- "H"
rownames(Alleles) <- rownames(Dia_freqDifs)


# save.image("x_fold_LD.Rdata")

#X-fold test - warning! I did not get around to parallelizing this script (yes I'm lazy) so each execution below can 
#take a couple of hours. 

chr1_LDH <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "1", LD = "H")
chr1_LDM <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "1", LD = "M")
chr1_LDL <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "1", LD = "L")
chr2_LDM <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "2", LD = "M")
chr2_LDL <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "2", LD = "L")
chr3_LDH <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "3", LD = "H")
chr3_LDM <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "3", LD = "M")
chr3_LDL <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "3", LD = "L")
chr4_LDH <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "4", LD = "H")
chr4_LDM <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "4", LD = "M")
chr4_LDL <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "4", LD = "L")
# chr5_LDH <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
#                     chr_exp_id = Alleles, chromo = "5", LD = "H")
# chr5_LDM <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
#                     chr_exp_id = Alleles, chromo = "5", LD = "M")
# chr5_LDL <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
#                     chr_exp_id = Alleles, chromo = "5", LD = "L")

#Separate run to lump all of the high LD groups from chromosome 2 into a single high LD group

ind <- grep("H"[1-9], Alleles$LDgr)
head(Alleles$LDgr[Alleles$LDgr == H[1-9]])
Alleles$LDgr[ind] <- "H"

chr2_LDH <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate, 
                    chr_exp_id = Alleles, chromo = "2", LD = "H")

###Rerun for the changes to chromosome 5 to account for potential sex determining region on chromosome 5. 

CDNDvSD <- read.table("newChr5_SDvCDND.txt", header=T)
Alleles$Dia <- CDNDvSD$freqDif

mer_Ids <- read.csv("/media/raglandlab/ExtraDrive4/DiaClass/clines/RAD_Library_Info_PostSequencing_FixedForPlateSwitch_MDids.csv", header=T)
female_Ids <- read.table("/media/raglandlab/ExtraDrive4/DiaClass/clines/ndsdcd.f.txt", header=F)
female_Ids <- as.character(female_Ids$V1)
SD_f_ids <- female_Ids[grep("sd", female_Ids)]
DiaND_f_ids <- female_Ids[grep("sd", female_Ids, invert = T)]
SD_f_ids <- as.character(mer_Ids$EcoRIAdapterName[mer_Ids$MD.ids %in% SD_f_ids])
DiaND_f_ids <- as.character(mer_Ids$EcoRIAdapterName[mer_Ids$MD.ids %in% DiaND_f_ids])
SD <- SD[,colnames(SD) %in% SD_f_ids]
CDND <- CDND[,colnames(CDND) %in% DiaND_f_ids]

chr5_LDH <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "5", LD = "H")
chr5_LDM <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "5", LD = "M")
chr5_LDL <- x_fold2(Dia1 = SD, Dia2 = CDND, Ecl1 = appleEarly, Ecl2 = appleLate, Ecl3 = HawEarly, Ecl4 = HawLate,
                    chr_exp_id = Alleles, chromo = "5", LD = "L")



#save.image("x_fold_chr2HLD.Rdata")

# load("x_fold_LD.Rdata")
# load("x_fold_chr2HLD.Rdata")

x_fold_est <- c(chr1_LDH[[1]][2], chr1_LDM[[1]][2], chr1_LDL[[1]][2], 
            chr2_LDH[[1]][2], chr2_LDM[[1]][2], chr2_LDL[[1]][2], 
            chr3_LDH[[1]][2], chr3_LDM[[1]][2], chr3_LDL[[1]][2], 
            chr4_LDH[[1]][2], chr4_LDM[[1]][2], chr4_LDL[[1]][2],
            chr5_LDH[[1]][2], chr5_LDM[[1]][2], chr5_LDL[[1]][2])

x_fold_pval <- c(chr1_LDH[[1]][1], chr1_LDM[[1]][1], chr1_LDL[[1]][1], 
                chr1_LDH[[1]][1], chr1_LDM[[1]][1], chr1_LDL[[1]][1], 
                chr3_LDH[[1]][1], chr3_LDM[[1]][1], chr3_LDL[[1]][1], 
                chr4_LDH[[1]][1], chr4_LDM[[1]][1], chr4_LDL[[1]][1],
                chr5_LDH[[1]][1], chr5_LDM[[1]][1], chr5_LDL[[1]][1])

#Print the raw data from the xfold tests 

chrs <- c(rep("1",3), rep("2", 3), rep("3",3), rep("4", 3), rep("5",3))
LD <-c(rep(c("high", "int", "low"), 5))
xfold_df <- as.data.frame(cbind(x_fold_est, x_fold_pval, chrs, LD))
colnames(xfold_df) <- c("x_fold", "p_val", "chromosomes", "LD")
write.table(xfold_df, file = "./results/x_fold_estimates.txt", quote = F, sep = "\t", row.names = F, col.names = T)

xfold_df$x_fold <- as.numeric(as.character(xfold_df$x_fold))
val <- c("*", "*", "", "*", "*", "","*", "*", "*", "", "", "", "*", "", "" )
library(ggplot2)
library(RColorBrewer)
library(viridis)
mypalette <- viridis_pal(option = "B")(1000)
p<- ggplot(data=xfold_df, aes(x=chromosomes, y=x_fold, fill=LD))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic() + xlab("chromosomes") +
  ylab("x-fold enrichments")
pdf("./results/x_fold_bar_new.pdf")
p + scale_fill_manual(values = c(mypalette[250], mypalette[500], mypalette[750]))+
  theme(axis.text = element_text(size = 14.0), axis.title = element_text(size=14.0))+
  geom_text(aes(label = val), vjust = -0.5, position = position_dodge(0.9), size=16)+
  geom_hline(yintercept=1)
dev.off()
