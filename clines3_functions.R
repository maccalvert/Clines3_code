##### Functions allele frequency difference permutation tests and allele frequeny difference correlation permutation tests. 

#Split the phred scaled genotype liklihoods and compute a mean genotype. 
splitPls <- function(x){
  a <-as.numeric(unlist(strsplit(x,","))) #splits into vectro of three characters, then converts to numeric
  b <- 10^(-a/10)
  c <- 1/(sum(b))
  d <- c*b
  e <- sum(d*c(0,1,2))
  return(e)
}

# empirical alternate allele frequency estimate. 
freqEst <- function(f){
  # If a locus has all NAs, then we assign it a zero - this will occasially happen during iteration sof the permutation test
  if(sum(is.na(f)) == length(f)) { 
    nf = 1
    f = 0
  } else {
    #remove NAs
    f <- f[!is.na(f)]
    nf = length(f)
  }
  #sum mean genotype and divide by number of individauls with data. This is multiplied by 2 in order to center the estimate on the alterate allele.
  freq<-(sum(f)/(2*nf))
  return(freq)
}

#Empirical allele frequency estimate, uses the frequency estimate function above
freqDif<-function(f1,f2) {
  dif<-freqEst(f1) - freqEst(f2)
  return(dif)}

#Function foe empirical allele frequency estimate between diapause termination phenotpes (early/late) over both apple and haw. 
#frequency estimates of early/late bulks within host races and then averaged across host races.
freqDif_2<-function(f1, f2, f3, f4) {
  dif <- mean(c(freqEst(f1), freqEst(f3))) - mean(c(freqEst(f2), freqEst(f4)))
  return(dif)
}

#Function for calculation of p-value using distribution of permutations and an empirical estimate. 
permPval<-function(est,vec) {
  prob<-ecdf(vec)(est)
  if (prob <= 0.5) {p=prob} else {p=1-prob}
  p=p*2 #accounts for 2-tailed search
  return(p)
}


#function for randomly sampling mean genotypes at a given locus to perform the permutation test for allele frequency differences.
#Prints p-value of the permutation test. 

permTest <- function(dia1, dia2, nit) {
  a <- length(dia1)
  b <- length(dia2)
  difs <- vector(length=nit)
  for(i in 1:nit){
    all <- sample(c(dia1,dia2))
    newdia <- all[1:a]
    newdia2 <- all[(a+1):(a+b)]
    difs[i] <- freqDif(newdia, newdia2)
  }
  pointEst<-freqDif(dia1, dia2)
  pval<-permPval(pointEst,difs)
  return(pval)
}

# Permutation test for significance of allele frequency differences between early and late diapause termination groups. To follow the 
# methods of Ragland et al. 2017, the early bulks and late bulks of each host race are average together to produce an average early and average late 
# allele frequency difference that are then subtracted from each other to form a average allele frequency difference between early and late 
# diapause termination.
#Dia 1 and Dia 3 must be the early bulks, Dia 2 and Dia 4 must be the late bulks. 
permTest_2 <- function(dia1, dia2, dia3, dia4, nit) {
  a <- length(dia1)
  b <- length(dia2)
  c <- length(dia3)
  d <- length(dia4)
  difs <- vector(length=nit)
  for(i in 1:nit){
    all <- sample(c(dia1,dia2))
    newdia <- all[1:a]
    newdia2 <- all[(a+1):(a+b)]
    newdia3 <- all[(a+b+1):(a+b+c)]
    newdia4 <- all[(a+b+c+1):(a+b+c+d)]
    early <- mean(c(freqEst(newdia), freqEst(newdia3)))
    late <- mean(c(freqEst(newdia2), freqEst(newdia4)))
    difs[i] <- early - late
  }
  pointEst <- freqDif_2(dia1, dia2, dia3, dia4)
  #pointEst<- mean(c(freqEst(dia1), freqEst(dia3))) - mean(c(freqEst(dia2), freqEst(dia4)))
  pval<-permPval(pointEst,difs)
  return(pval)
}

#Permutation test functions for test of significant snp enrichment 
permTest_perc <- function(dia1, dia2, nit) {
  a <- length(dia1)
  b <- length(dia2)
  difs <- vector(length=nit)
  for(i in 1:nit){
    all <- sample(c(dia1,dia2))
    newdia <- all[1:a]
    newdia2 <- all[(a+1):(a+b)]
    difs[i] <- freqDif(newdia, newdia2)
  }
  return(difs)
}


#Permutation test for over representation of signficant SNPs on specific genomic regions 
#emp_AF is a table of allele frequencies, it must have a column called "p_val" with p values
# permPercent <- function(dia1, dia2, nit, emp_AF){
#   pval_mat <- matrix(nrow=nit, ncol=nrow(dia1))
#   for(i in 1:nrow(dia1)){
#     a <- length(dia1[i,])
#     b <- length(dia2[i,])
#     difs <- vector(length=nit)
#     for(j in 1:nit){
#       all <- unlist(sample(c(dia1[i,],dia2[i,])))
#       newdia <- all[1:a]
#       newdia2 <- all[(a+1):(a+b)]
#       difs[j] <- freqDif(newdia, newdia2)
#     }
#     pval_mat[,i] <- difs
#   }
#   sums <- vector(length=nrow(pval_mat))
#   for(i in 1:nrow(pval_mat)){
#     for(j in 1:ncol(pval_mat)){
#       pval_mat[i,j] <- permPval(pval_mat[i,j], pval_mat[,j])
#     }
#   }
#   for(i in 1:nrow(pval_mat)){
#     sums[i] <- sum(pval_mat[i,] < 0.05)/length(pval_mat[i,])
#   }
#   emp_sum <- sum(emp_AF$p_val < 0.05)/length(emp_AF$p_val)
#   emp_sum_pval <- permPval(emp_sum, sums)
#   sum_r <- c(emp_sum, emp_sum_pval)
#   names(sum_r) <- c("emp_perc", "p_val")
#   return(sum_r)
# }


#Returns a permuted allele frequency difference for use in the correlation permutation test in the 
#function below. 
permForCor <- function(dia1, dia2) {
  a <- length(dia1)
  b <- length(dia2)
  all <- sample(c(dia1,dia2))
  newdia <- all[1:a]
  newdia2 <- all[(a+1):(a+b)]
  dif <- freqDif(newdia, newdia2)
  return(dif)
}


#permutation function for correlations.
permCor <- function(dia1, dia2, ecl1, ecl2, nit){
  vec <- vector(length=nit)
  for(i in 1:nit){
    rep <- vector(length=nrow(dia1))
    for(j in 1:nrow(dia1)){
      rep[j] <- permForCor(unlist(dia1[j,]), unlist(dia2[j,]))
    }
    rep_2 <- vector(length=nrow(ecl1))
    for(k in 1:nrow(ecl1)){
      rep_2[k]<- permForCor(unlist(ecl1[k,]), unlist(ecl2[k,]))
    }
    a<-cor.test(rep, rep_2, method = "spearman")
    vec[i] <- a$estimate
  }
  return(vec)
}

#This function is meant to deal with the treatment of loci on chromosome 5, for which we only evaluate the difference. This matters for 
#later analysis of chromosome and LD groups that will essentially involve differnet individuals for different parts of the analysis 
# among female 
permCorForLDs <- function(dia1, dia2, ecl1, ecl2, chr5_1, chr5_2, nit){
  vec <- vector(length=nit)
  for(i in 1:nit){
    rep <- vector(length=nrow(dia1))
    for(j in 1:nrow(dia1)){
      rep[j] <- permForCor(unlist(dia1[j,]), unlist(dia2[j,]))
    }
    names(rep) <- rownames(dia1)
    rep_2 <- vector(length=nrow(ecl1))
    for(k in 1:nrow(ecl1)){
      rep_2[k]<- permForCor(unlist(ecl1[k,]), unlist(ecl2[k,]))
    }
    names(rep_2) <- rownames(ecl1)
    for(l in (nrow(dia1)+1):(nrow(dia1)+nrow(chr5_1))){
      rep[l] <- permForCor(unlist(chr5_1[(l-nrow(dia1)),]), unlist(chr5_2[(l-nrow(dia1)),]))
    }
    names(rep)[(nrow(dia1)+1):(nrow(dia1)+nrow(chr5_1))] <- rownames(chr5_1)
    rep <- rep[order(match(names(rep), names(rep_2)))]
    if(identical(names(rep), names(rep_2)) == T){
      a<-cor.test(rep, rep_2, method = "spearman")
      vec[i] <- a$estimate
    }else{
      stop("rownames of combined permuted allele frequeny tables do not match, check input files")
    }
  }
  return(vec)
}

#For the correlation permutation tests that involve the emergence timing data. 
permForCor_2 <- function(dia1, dia2, dia3, dia4) {
  a <- length(dia1)
  b <- length(dia2)
  c <- length(dia3)
  d <- length(dia4)
  all <- sample(c(dia1,dia2,dia3,dia4))
  newdia <- all[1:a]
  newdia2 <- all[(a+1):(a+b)]
  newdia3 <- all[(a+b):(a+b+c)]
  newdia4 <- all[(a+b+c+1):(a+b+c+d)]
  early <- mean(c(freqEst(newdia), freqEst(newdia3)))
  late <- mean(c(freqEst(newdia2), freqEst(newdia4)))
  dif <- early - late
  return(dif)
}

permReg_2 <- function(dia1, dia2, ecl1, ecl2, ecl3, ecl4, nit){
  vec <- vector(length=nit)
  for(i in 1:nit){
    rep <- vector(length=nrow(dia1))
    for(j in 1:nrow(dia1)){
      rep[j] <- permForCor(unlist(dia1[j,]), unlist(dia2[j,]))
    }
    rep_2 <- vector(length=nrow(ecl1))
    for(k in 1:nrow(ecl1)){
      rep_2[k]<- permForCor_2(unlist(ecl1[k,]), unlist(ecl2[k,]), unlist(ecl3[k,]), unlist(ecl4[k,]))
    }
    a<-cor.test(rep, rep_2, method = "spearman")
    vec[i] <- a$estimate
  }
  return(vec)
}

permCorForLDs_2 <- function(dia1, dia2, ecl1, ecl2, ecl3, ecl4, chr5_1, chr5_2, nit){
  vec <- vector(length=nit)
  for(i in 1:nit){
    rep <- vector(length=nrow(dia1))
    for(j in 1:nrow(dia1)){
      rep[j] <- permForCor(unlist(dia1[j,]), unlist(dia2[j,]))
    }
    names(rep) <- rownames(dia1)
    rep_2 <- vector(length=nrow(ecl1))
    for(k in 1:nrow(ecl1)){
      rep_2[k]<- permForCor_2(unlist(ecl1[k,]), unlist(ecl2[k,]), unlist(ecl3[k,]), unlist(ecl4[k,]))
    }
    names(rep_2) <- rownames(ecl1)
    for(l in (nrow(dia1)+1):(nrow(dia1)+nrow(chr5_1))){
      rep[l] <- permForCor(unlist(chr5_1[(l-nrow(dia1)),]), unlist(chr5_2[(l-nrow(dia1)),]))
    }
    names(rep)[(nrow(dia1)+1):(nrow(dia1)+nrow(chr5_1))] <- rownames(chr5_1)
    rep <- rep[order(match(names(rep), names(rep_2)))]
    if(identical(names(rep), names(rep_2)) == T){
      a<-cor.test(rep, rep_2, method = "spearman")
      vec[i] <- a$estimate
    }else{
      stop("rownames of combined permuted allele frequeny tables do not match, check input files")
    }
  }
  return(vec)
}

#Asignments are strange because of the way loci are parraleled, essentially we're asking if 
#Ecl and diapause alleles show up in the same quadrat of the cartesian plane for often than 
# would be expected by chance 

permTest_XF <- function(dia1, dia2) {
  a <- length(dia1)
  b <- length(dia2)
  all <- sample(c(dia1,dia2))
  newdia <- all[1:a]
  newdia2 <- all[(a+1):(a+b)]
  difs <- freqDif(newdia, newdia2)
  return(difs)
}

permTest_2_XF <- function(dia1, dia2, dia3, dia4) {
  a <- length(dia1)
  b <- length(dia2)
  c <- length(dia3)
  d <- length(dia4)
  all <- sample(c(dia1,dia2))
  newdia <- all[1:a]
  newdia2 <- all[(a+1):(a+b)]
  newdia3 <- all[(a+b+1):(a+b+c)]
  newdia4 <- all[(a+b+c+1):(a+b+c+d)]
  early <- mean(c(freqEst(newdia), freqEst(newdia3)))
  late <- mean(c(freqEst(newdia2), freqEst(newdia4)))
  difs <- early - late
  return(difs)
}


x_fold2 <- function(Dia1, Dia2, Ecl1, Ecl2, Ecl3, Ecl4, chr_exp_id, chromo, LD){
  #calculates the observed value for each locus for the desired chromosome / LD group
  chr_exp <- as.data.frame(chr_exp_id)
  chr_exp <- chr_exp[chr_exp$chr == chromo,]
  chr_exp <- chr_exp[complete.cases(chr_exp),]
  chr_exp <- chr_exp[chr_exp$LDgr==LD,]
  dia <- chr_exp$Dia
  ecl <- chr_exp$Ecl
  dia[dia > 0] = 1
  dia[dia < 0] = -1
  ecl[ecl > 0] = 1
  ecl[ecl < 0] = -1
  chrsigns = dia * ecl
  obsr2 = sum((dia*ecl)==-1)
  #Now we need to subset each of the genotype files by the chromosome / LD group 
  ind <- rownames(Dia1) %in% rownames(chr_exp)
  Dia1 <- Dia1[ind,]
  Dia2 <- Dia2[ind,]
  Ecl1 <- Ecl1[ind,]
  Ecl2 <- Ecl2[ind,]
  Ecl3 <- Ecl3[ind,]
  Ecl4 <- Ecl4[ind,]
  #Now calculate the Null distribution of permuted genotypes for this chromosome 
  N <- length(1:nrow(chr_exp))
  null = rep(NA, 10000)
  int <- vector(length=nrow(Dia1))
  ecl <- vector(length=nrow(Ecl1))
  #This is slow, we need to find a way to paralellize it
  # The way it works now is by calculating permuted allele frequency differences for each SNP 
  # diapause intensity and diapause termination experiments. Then Extracts the "sign" (positive or negative) of those permutations 
  # Then sums the number of loci in each experiment that displays sign concordance - i.e. the permutated allele frequency difference 
  # has the same sign in both experiments. This is all done 10,000 times. The p value us calculated by determining what proportion 
  # of iterations have a value greater than the empirical observation. The x fold statistic is caluclated as the empirical estimate 
  # over the mean of the null distribution. 
  for (i in 1:10000){
    for (j in 1:nrow(Dia1)){
      int[j] <- permTest_XF(dia1 = unlist(Dia1[j,]), dia2 = unlist(Dia2[j,]))
      ecl[j] <- permTest_2_XF(dia1 = unlist(Ecl1[j,]), dia2 = unlist(Ecl2[j,]), dia3 = unlist(Ecl3[j,]), dia4 = unlist(Ecl4[j,]))
    }
    int[int > 0] = 1
    int[int < 0] = -1
    ecl[ecl > 0] = 1
    ecl[ecl < 0] = -1
    null_sign <- int*ecl
    null[i] = sum(null_sign == -1)
  }
  p2 = mean(null > obsr2)
  xfold2 = obsr2 / mean(null)
  resultsdat = cbind(p2, xfold2, obsr2)
  return(list(resultsdat, null))
}

#This function polarizes loci according to an index that says the is (T) or is not (F) to be polarized
pol_func <- function(vec, pol){
  for(i in 1:length(vec)){
    if(pol[i] == T){
      vec[i] <- 1-vec[i]
    }else if(pol[i] == F){
      vec[i] <- vec[i]}}
  return(vec)}

#Calculation of the 95% confidence interval for mean allele frequency estimates derived from sets of genotyped loci. 
CI_95 <- function(vec, n){
  se <- sd(vec)/sqrt(n)
  Q <- qnorm(1- 0.05/2)
  interval <- c(mean(vec)-Q*se, mean(vec)+Q*se )
  return(interval)
}


#############################################
#Earlier iterations of the functions above.##
#############################################

# permTest <- function(dia1, dia2) {
#   #determine number of individuals in each group
#   a <- length(dia1)
#   b <- length(dia2)
#   #Randomly shuffle genotypes
#   all <- sample(c(dia1,dia2))
#   #assign to new groups with same size as original groups
#   newdia <- all[1:a]
#   newdia2 <- all[(a+1):(a+b)]
#   #Estimate the new allele frequency difference
#   dif <- freqDif(newdia, newdia2)
#   return(dif)
# }
# 
# #Produces a vector of permuted allele frequences that can then be fed into the pvalue function below
# #permutation tests retains unbalanced sample sizes for both groups. Necessary for following
# #assumption that values in each group are independent and identicall distributed and for creating 
# #a null that is refleftive of the information in the group means. 
# 
# permute <- function(dia1, dia2, nit) {
#   a <- length(dia1)
#   b <- length(dia2)
#   difs <- vector(length=nit)
#   for(i in 1:nit){
#     all <- sample(c(dia1,dia2))
#     newdia <- all[1:a]
#     newdia2 <- all[(a+1):(a+b)]
#     difs[i] <- freqDif(newdia, newdia2)
#   }
#   return(difs)
# }
# 

