#################
### NEP.calc()
#################

#Description
#This funciton calculations non-exclusion power for any given set of loci with any given set of alleles
#originally from Jamieson, A. 1997

#Input Parameters:
#tmp - a data.frame with only two columns, the first is the name of the loci and the second the frequencies
#      of the alleles at that locus - in long form

nep.calc <- function(tmp, opt=1){
  loci <- unique(tmp[,1])

  i <- NULL
  OUT <- NULL
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- -2*sum(tmp1^2)
    p2 <- sum(tmp1^3)
    p3 <- 2*sum(tmp1^4)
    p4 <- -3*sum(tmp1^5)
    p5 <- -2*(sum(tmp1^2)^2)
    p6 <- 3*sum(tmp1^2)*sum(tmp1^3)
    P <- 1+(sum(p1,p2,p3,p4,p5,p6))
    OUT <- c(OUT,P)
  }

  i<-1
  i <- NULL
  OUT1 <- NULL
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- -4*sum(tmp1^2)
    p2 <- 2*(sum(tmp1^2)^2)
    p3 <- 4*sum(tmp1^3)
    p4 <- -3*sum(tmp1^4)
    P <- 1+(sum(p1,p2,p3,p4))
    OUT1 <- c(OUT1,P)
  }

  i<-1
  i <- NULL
  OUT2 <- NULL
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- 4*sum(tmp1^4)
    p2 <- -4*sum(tmp1^5)
    p3 <- -3*sum(tmp1^6)
    p4 <- -8*(sum(tmp1^2)^2)
    p5 <- 8*sum(tmp1^2)*sum(tmp1^3)
    p6 <- 2*(sum(tmp1^3)^2)
    P <- 1+(sum(p1,p2,p3,p4,p5,p6))
    OUT2 <- c(OUT2,P)
  }

  #parparing for first for loops
  i <-  NULL
  PId <- NULL
  for(i in 1:length(loci)){

    #grabbing each locus individually and ordering some smallest to largest
    tmp1 <- tmp[tmp[,1] == loci[i],]
    tmp1 <- tmp1[order(tmp1[,2]),]

    #raising the allele frequencies to the fourth power
    p14 <-  sum(tmp1[,2]^4)
    p14 #makes test in excel

    #parparing for second for loops
    j <- NULL
    calc <- NULL

    for(j in 1:nrow(tmp1)){

      #take each allele frequencies one at a time
      tmp1.f <- tmp1[j,2]
      tmp1.f

      #and taking all allele frequencies greater than that and (2*Pi*Pj)^2 and summing
      tmp2 <- tmp1[-1:-j,2]
      calc1 <- sum((2*tmp1.f*tmp2)^2)

      #adding to calc
      calc <- c(calc,calc1)
    }

    #summing all calcs
    calc <- sum(calc)

    #adding p14 and calc
    PId1 <- p14 +calc

    #adding to PId
    PId <- c(PId,PId1)
  }

  #parparing for first for loops
  i <-  NULL
  PrSib <- NULL
  for(i in 1:length(loci)){

    #grabbing each locus individually and ordering some smallest to largest
    tmp1 <- tmp[tmp[,1] == loci[i],]
    tmp1 <- tmp1[order(tmp1[,2]),]
    tmp1$p2 <- tmp1[,2]^2
    tmp1$p4 <- tmp1[,2]^4
    p1 <- 0.25
    p2 <- 0.5*sum(tmp1$p2)
    p3 <- 0.5*(sum(tmp1$p2)^2)
    p4 <- -0.25*(sum(tmp1$p4))
    c(p1,p2,p3,p4)
    ps1 <- sum(p1,p2,p3,p4)
    PrSib <- c(PrSib,ps1)
  }

  #getting the power calcs accross all loci and adding it to the individual locus calcs
  Locus <- c(loci,"Combined")

  NE_2P <- c(OUT , prod(1-OUT))
  NE_1P <- c(OUT1 , prod(1-OUT1))
  NE_PP <- c(OUT2 , prod(1-OUT2))
  Pid <- c(PId,prod(PId))
  PrSib <- c(PrSib,prod(PrSib))
  Output <- data.frame(Locus,NE_1P,NE_2P,NE_PP, Pid,PrSib, stringsAsFactors=F)

  if(opt==1){
    return(Output)
  }

  if(opt==2){
    Output <- Output[Output$Locus == "Combined",]
    return(Output)
  }
}
