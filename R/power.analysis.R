########################
### power.analysis() ###
########################

#Description
#This fuction is a modified version of Mark Christie's script to calculate the exclusion power within a dataset
# It will calculates the expected number of false parent offspring pairs (POPs) within the dataset and
# it also calculates the exclusion power per locus


#Input Parameters:
#adults - this is a data frame with parents ids and genotypes
#offspring - this is a data frame with offspring ids and genotypes
#loci_number - define as the number of loci you want to analyze
#             NOTE: from what I understand it randomly selects which ones to analyze
#opt - has three versions
#       opt 1 - writes both Exp. Number of false POPs and Locus.Power to file
#       opt 2 - returns Exp. Number of false POPs
#       opt 3 - reutnrs Locus.Power

power.analysis <- function(adults,offspring, loci_number=0,opt=1){

  if(length(loci_number)==1){
    if(loci_number == 0){print("Need to define the number of loci to analyze before I can continue..."); stop}
    multi <- T
  } else {
    multi <- F
  }

  if(multi==F){

    #a long winded verson of getting the adults and offspring gts without their sample names
    Adults1 <- adults
    a <- ncol(Adults1)
    lociset <- seq(from= 2,to= a,by= 2)
    loci <- as.numeric(loci_number)
    s1 <- sample(lociset,loci,replace= FALSE)
    s2 <- sort(c(s1,s1+1))
    Adults <- Adults1[,s2]
    Offspring1 <- offspring
    Offspring <- Offspring1[,s2]

    head(Offspring)
    head(Adults)
    a <- 1
    #Begin calcualtion of Pr(Z) and Pr(delta)=======================================#
    massive <- function(massive){

      #making a table with unique alleles, their counts, and frequencies for parents
      AAT40 <- c(Adults[,i],Adults[,(i+1)])
      locus <- as.data.frame(table(AAT40))
      Frequency <- locus[,2]/sum(locus[,2])
      Data <- cbind(locus, Frequency)
      Data

      #removing missing data from the table
      if (Data[1,1]==0) {

        Data <- Data[-1,]

      } else {

        Data <- Data
      } #end of missing data remover


      #calculating the expected number of each homozygote allele
      n1 <- ((length(AAT40))/2)
      A1 <- (Data[,3]*(2*n1))
      A2 <- (Data[,3]^2)*n1
      AA <- A1-A2
      AAA <- cbind(Data,AA)
      AAA

      #making a table with unique alleles, their counts, and frequencies for offspring
      AAT402 <- c(Offspring[,i],Offspring[,(i+1)])
      locus2 <- as.data.frame(table(AAT402))
      Frequency2 <- locus2[,2]/sum(locus2[,2])
      Data2 <- cbind(locus2, Frequency2)
      Data2

      #removing missing data again
      if (Data2[1,1]==0){

        Data2=Data2[-1,]

      } else {

        Data2=Data2
      } #end of missing data remover

      #calculating the expected number of each homozygote allele
      n2 <- ((length(AAT402))/2)
      B1 <- (Data2[,3]*(2*n2))
      B2 <- (Data2[,3]^2)*n2
      BB <- B1-B2
      BBB <- cbind(Data2,BB)
      BBB

      #getting the alleles in the parents that are found in the kids
      AAA1 <- match(BBB[,1],AAA[,1])
      g <- which(AAA1>0)
      g1 <- AAA1[g]
      AAA1 <- AAA[g1,]
      AAA1

      #doing the same thing for the kids
      BBB1 <- match(AAA[,1],BBB[,1])
      j <- which(BBB1>0)
      j1 <- BBB1[j]
      BBB1 <- BBB[j1,]

      #multiplying the number of expected from parents and kids and summing
      AB <- AAA1[,4]*BBB1[,4]
      AB <- sum(AB)
      AB

      #getting all possible combinations of alleles and removing homozygotes for the parents
      y <- AAA1[,1]
      Agenotypes <- expand.grid(y,y)
      remove <- which(Agenotypes[,1]==Agenotypes[,2])
      Agenotypes <- Agenotypes[-remove,]
      Agenotypes

      #doing the same for the kids
      z <- BBB1[,1]
      Bgenotypes <- expand.grid(z,z)
      remove <- which(Bgenotypes[,1]==Bgenotypes[,2])
      Bgenotypes <- Bgenotypes[-remove,]
      Bgenotypes

      #getting the expected number of parents that would have each heterozygous genotype
      one <- match(Agenotypes[,1],AAA1[,1])
      two <- match(Agenotypes[,2],AAA1[,1])
      one <- AAA1[one,3]
      two <- AAA1[two,3]
      Agfreq <- one*two*n1*2
      Agfreq

      #getting the expected number of kids that would have each heterozygous genotype
      oneb <- match(Bgenotypes[,1],BBB1[,1])
      twob <- match(Bgenotypes[,2],BBB1[,1])
      oneb <- BBB1[oneb,3]
      twob <- BBB1[twob,3]
      Ogfreq <- oneb*twob*n2*2
      Ogfreq

      #combining the two and calculating the probability of being false at that locus
      Gfreqs <- Agfreq*Ogfreq
      Gfreqs <- floor(Gfreqs)
      Gfreq <- sum(Gfreqs)/2
      PrB <- (AB-Gfreq)/(n1*n2)
      write.table(PrB,file="Output_Per_Locus_Exclusion_Probabilities.txt",row.names=FALSE,col.names=F,sep="\t",append=T)

    } # end of massive function

    a <- ncol(Adults)
    loci.locs <- seq(from=1,to=ncol(Adults)-1, by=2)
    loci.locs

    #applying massive function to all loci
    for(i in loci.locs) {
      lapply(i,massive)
    } # end of for loop calculating exclusion probs

    #reading the file with each locus' exclusion prob in
    PrBs <- read.table("Output_Per_Locus_Exclusion_Probabilities.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

    #getting the number of parents and offspring
    nn1 <- length(Adults[,1])
    nn2 <- length(Offspring[,1])

    #getting the product of all the independent probs
    Pr_delta <- prod(PrBs[,1])

    #getting exclusion power per locus
    Locus.Power <- data.frame(colnames(Offspring)[loci.locs],PrBs,stringsAsFactors=F)
    colnames(Locus.Power) <- c("Loci","Exc.Power")
    Locus.Power

    #making table with information in it
    Expected_Number_of_False_Pairs <- Pr_delta*(nn1*nn2)
    pvalue1 <- data.frame(loci_number,Expected_Number_of_False_Pairs,Pr_delta,stringsAsFactors=F)

    #removing the saved file
    unlink("Output_Per_Locus_Exclusion_Probabilities.txt")

  } else {

    pvalue1 <- NULL
    for(i in 1:length(loci_number)){

      #a long winded verson of getting the adults and offspring gts without their sample names
      Adults1 <- adults
      a <- ncol(Adults1)
      lociset <- seq(from= 2,to= a,by= 2)
      loci <- as.numeric(loci_number[i])
      s1 <- sample(lociset,loci,replace= FALSE)
      s2 <- sort(c(s1,s1+1))
      Adults <- Adults1[,s2]
      Offspring1 <- offspring
      Offspring <- Offspring1[,s2]

      head(Offspring)
      head(Adults)
      a <- 1
      #Begin calcualtion of Pr(Z) and Pr(delta)=======================================#
      massive <- function(massive){

        #making a table with unique alleles, their counts, and frequencies for parents
        AAT40 <- c(Adults[,i],Adults[,(i+1)])
        locus <- as.data.frame(table(AAT40))
        Frequency <- locus[,2]/sum(locus[,2])
        Data <- cbind(locus, Frequency)
        Data

        #removing missing data from the table
        if (Data[1,1]==0) {

          Data <- Data[-1,]

        } else {

          Data <- Data
        } #end of missing data remover


        #calculating the expected number of each homozygote allele
        n1 <- ((length(AAT40))/2)
        A1 <- (Data[,3]*(2*n1))
        A2 <- (Data[,3]^2)*n1
        AA <- A1-A2
        AAA <- cbind(Data,AA)
        AAA

        #making a table with unique alleles, their counts, and frequencies for offspring
        AAT402 <- c(Offspring[,i],Offspring[,(i+1)])
        locus2 <- as.data.frame(table(AAT402))
        Frequency2 <- locus2[,2]/sum(locus2[,2])
        Data2 <- cbind(locus2, Frequency2)
        Data2

        #removing missing data again
        if (Data2[1,1]==0){

          Data2=Data2[-1,]

        } else {

          Data2=Data2
        } #end of missing data remover

        #calculating the expected number of each homozygote allele
        n2 <- ((length(AAT402))/2)
        B1 <- (Data2[,3]*(2*n2))
        B2 <- (Data2[,3]^2)*n2
        BB <- B1-B2
        BBB <- cbind(Data2,BB)
        BBB

        #getting the alleles in the parents that are found in the kids
        AAA1 <- match(BBB[,1],AAA[,1])
        g <- which(AAA1>0)
        g1 <- AAA1[g]
        AAA1 <- AAA[g1,]
        AAA1

        #doing the same thing for the kids
        BBB1 <- match(AAA[,1],BBB[,1])
        j <- which(BBB1>0)
        j1 <- BBB1[j]
        BBB1 <- BBB[j1,]

        #multiplying the number of expected from parents and kids and summing
        AB <- AAA1[,4]*BBB1[,4]
        AB <- sum(AB)
        AB

        #getting all possible combinations of alleles and removing homozygotes for the parents
        y <- AAA1[,1]
        Agenotypes <- expand.grid(y,y)
        remove <- which(Agenotypes[,1]==Agenotypes[,2])
        Agenotypes <- Agenotypes[-remove,]
        Agenotypes

        #doing the same for the kids
        z <- BBB1[,1]
        Bgenotypes <- expand.grid(z,z)
        remove <- which(Bgenotypes[,1]==Bgenotypes[,2])
        Bgenotypes <- Bgenotypes[-remove,]
        Bgenotypes

        #getting the expected number of parents that would have each heterozygous genotype
        one <- match(Agenotypes[,1],AAA1[,1])
        two <- match(Agenotypes[,2],AAA1[,1])
        one <- AAA1[one,3]
        two <- AAA1[two,3]
        Agfreq <- one*two*n1*2
        Agfreq

        #getting the expected number of kids that would have each heterozygous genotype
        oneb <- match(Bgenotypes[,1],BBB1[,1])
        twob <- match(Bgenotypes[,2],BBB1[,1])
        oneb <- BBB1[oneb,3]
        twob <- BBB1[twob,3]
        Ogfreq <- oneb*twob*n2*2
        Ogfreq

        #combining the two and calculating the probability of being false at that locus
        Gfreqs <- Agfreq*Ogfreq
        Gfreqs <- floor(Gfreqs)
        Gfreq <- sum(Gfreqs)/2
        PrB <- (AB-Gfreq)/(n1*n2)
        write.table(PrB,file="Output_Per_Locus_Exclusion_Probabilities.txt",row.names=FALSE,col.names=F,sep="\t",append=T)

      } # end of massive function

      a <- ncol(Adults)
      loci.locs <- seq(from=1,to=ncol(Adults)-1, by=2)
      loci.locs

      #applying massive function to all loci
      for(i in loci.locs) {
        lapply(i,massive)
      } # end of for loop calculating exclusion probs

      #reading the file with each locus' exclusion prob in
      PrBs <- read.table("Output_Per_Locus_Exclusion_Probabilities.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

      #getting the number of parents and offspring
      nn1 <- length(Adults[,1])
      nn2 <- length(Offspring[,1])

      #getting the product of all the independent probs
      Pr_delta <- prod(PrBs[,1])

      #getting exclusion power per locus
      Locus.Power <- data.frame(colnames(Offspring)[loci.locs],PrBs,stringsAsFactors=F)
      colnames(Locus.Power) <- c("Loci","Exc.Power")
      Locus.Power

      #making table with information in it
      Expected_Number_of_False_Pairs <- Pr_delta*(nn1*nn2)
      pvalue1.1 <- data.frame(Expected_Number_of_False_Pairs,Pr_delta,stringsAsFactors=F)
      pvalue1 <- rbind(pvalue1, pvalue1.1)

      #removing the saved file
      unlink("Output_Per_Locus_Exclusion_Probabilities.txt")
    }
    pvalue1 <- cbind(loci_number,pvalue1)
  }

  #writing both to file
  if(opt == 1){

    write.table(Locus.Power,file="Locus power.txt",row.names=FALSE,col.names=T,sep="\t",append=T)
    write.table(pvalue1,file="Expected Number of False Pairs.txt",row.names=FALSE,col.names=T,sep="\t",append=T)

    print("Locus power.txt saved to current working directory")
    print("Expected Number of False Pairs.txt saved to current working directory")
  } #end of option 1

  #returning only the expected number of false pairs
  if(opt == 2){

    return(pvalue1)
  } #end of option 2

  #returning the exc. power per locus
  if(opt == 3){

    return(Locus.Power)
  } #end of option 3

} #end of power analysis function
