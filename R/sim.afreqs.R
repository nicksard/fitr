####################
### sim.afreqs() ###
####################

#Description
#This fuction is a modified version of Mark Christie's script to create simulate genotype allele frequencies.
#Note that there are high and low frequencies for each locus, and the alleles and their associated frequencies
#are exactly the same.

#Input Parameters:
#nlcoi - define the number of loci you want
#n.alleles - define the number of alleles you want


sim.afreqs <- function(nloci=10,n.alleles=10, simple=F){

  #making sure nloci and n.alleles are numeric
  nloci <- as.numeric(nloci)
  n.alleles <- as.numeric(n.alleles)

  #getting the alleles frequencies for each allele
  a <- c(1:n.alleles)
  b <- rev(a)
  c <- b/a
  d <- sum(c)

  #now making actual allele frequencies
  freqs <- c/d

  #setting the lowest allele to an arbitrary number
  lowest.allele <- 140

  #making afreqs

  #getting loci names
  loci <- paste0("Locus.",1:nloci)
  loci <- rep(loci, each=n.alleles)

  #for alternative output
  loci.1 <- rep(1:nloci, each=n.alleles)

  #getting the allele names
  alleles2 <- cbind(seq(lowest.allele,lowest.allele+n.alleles-1,1),freqs)
  alleles <- rep(alleles2[,1],nloci)


  #getting their frequencies
  freqs  <-  rep(alleles2[,2],nloci)

  #returning finished product
  if(simple == F){

    afreqs <- data.frame(loci,alleles,freqs, stringsAsFactors=F)
    afreqs
    return(afreqs)

  } else {

    afreqs <- data.frame(loci.1,freqs, stringsAsFactors=F)
    afreqs
    return(afreqs)
  }


} # end of sim.afreqs
