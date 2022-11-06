











#DataGeneration_TPD(Beta1=0.05,Sigma2=1,LD.r=0.75,P_g=0.2,P_z=0.3)


# N: phase 1 sample size, i.e. GWAS data
# Beta0: intercept
# Beta1: genetic effect
# Sigma2: variance of the error term 
# LD.r: linkage disequilibrium (r) between G and Z
# P_g: minor allele frequency for G, the causal SNP
# P_z: minor allele frequency for Z, the GWAS SNP
# tao: quantile value to define the stratification for the quantitative trait




#' This code a function that generates a sample Y, Z, G given some initial parameters.
#'
#'
#' @param Beta0  intercept  Default 2
#' @param Beta1  genetic effect   Default 0.5
#' @param Sigma2  variance of the error term   Default 1
#' @param N Phase 1 sample size, i.e. GWAS data (Default: 5000)
#' @param LD.r linkage disequilibrium (r) between G and Z  (Default: 0.75)
#' @param P_g  minor allele frequency for G, the causal SNP
#' @param P_z  minor allele frequency for Z, the GWAS SNP
#' @param tao quantile value to define the stratification for the quantitative trait (default: 2/5)
#' 
#' @return  A dataframe with complete data Y, G, Z, S, where G and Z come from the same haplotype determined by P_g, P_z and LD.r; Y is generated from Y = Beta0 + Beta1 x G, S is a 3 level variable determined by Y, Beta1, Sigma2 and tao. the function iterates across generated datasets until Z and Y are associated at a suggestive genome wide threshold of p<=1e-05.
#' 
#' 
#' @examples
#'
#' data = DataGeneration_TPD()
#'
#' @export
DataGeneration_TPD <- function(Beta0 = 2, Beta1 = 0.5, Sigma2 = 1, N = 5000, LD.r = 0.75, P_g = 0.2, P_z = 0.3, tao = 2/5)  {
  
  
  cat("Beta0=",Beta0, " Beta1=", Beta1, "\n");
  
  
  # P_g<-0.2; P_z<-0.3 #MAF of G and Z respectively
  # tao = quantile to define the strata in Y
  
  P_G <- 1 - P_g; P_Z <- 1 - P_z ;
  
  Rsq <- LD.r^2
  LD <- LD.r * sqrt(P_G*P_g*P_Z*P_z)
  # Dpr <- LD /( (LD<0)*min(P_G*P_Z,P_g*P_z) + (LD>=0)*min(P_z*P_Z,P_G*P_g) )
  # Dpr <- LD /( (LD<0)*min(P_G*P_Z,P_g*P_z) + (LD>=0)*min(P_G*P_z,P_g*P_Z) ) ## corrected on Nov 2018 following Wikipedia and https://doi.org/10.1186/1471-2105-8-428
  Dpr <- LD /( (LD<0)*max(-P_G*P_Z,-P_g*P_z) + (LD>=0)*min(P_G*P_z,P_g*P_Z) ) ## corrected on Apr 2019 following Wikipedia and https://doi.org/10.1186/1471-2105-8-428
  
  #Values
  #cat(Rsq, LD, Dpr,"\n")
  
  ## Calculate haplotype frequencies (two SNPs)
  h.freqs <- rep(0, 4)
  h.freqs[1] <- LD + (1-P_g)*(1-P_z)
  h.freqs[2] <- 1- P_g - h.freqs[1]
  h.freqs[3] <- 1 - P_z - h.freqs[1]
  h.freqs[4] <- P_g - h.freqs[3]
  
  pval<-1
  k<-0
  while( pval>=1e-05 ){
    h1 <- sample(1:4, size=N, prob=h.freqs, replace=TRUE)
    h2 <- sample(1:4, size=N, prob=h.freqs, replace=TRUE)
    G <- as.numeric(h1==3 | h1==4) + as.numeric(h2==3 | h2==4)
    Z <- as.numeric(h1==2 | h1==4) + as.numeric(h2==2 | h2==4)
    #Simulate the data with G the seq SNP and Z the tag SNP
    Y <- (Beta0 + G*Beta1 + rnorm(N,sd=sqrt(Sigma2)))
    Yst <- sapply(Y,function(x){
      if( x<qnorm(tao,mean=Beta0,sd=sqrt(Sigma2)) ){return(1)
      }else if( x>qnorm(1-tao,mean=Beta0,sd=sqrt(Sigma2)) ){return(3)
      }else return(2)})
    lmfit.comp <- lm(Y~Z)   
    pval <- 1-pchisq(as.numeric(coef(lmfit.comp)[2]^2/diag(vcov(lmfit.comp))[2]),1)
    k <- k+1
    if( pval<1e-05 ){       
      G0 <- sample(0:2, size=N, prob=c(P_G^2,2*P_G*P_g,P_g^2), replace=TRUE)
      dat_sim_it <- data.frame(wait_it=k,Y=Y,G1=G,Z=Z,G0=G0,S=Yst)
    }
  }
  return(dat_sim_it)
}



