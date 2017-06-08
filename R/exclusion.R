#################
### exclusion()
#################

#Description
#This funciton takes a data.frame with parents and a data.frame with offspring and identify all parent offspring
#relationships with up to the number of mismatches you describe and it has two alternative options for output
#take from SOLOMON

#Input Parameters:
#adults - contains the sample names and genotypes of the parents
#offspring - contains the sample names and genotypes of the offspring
#mismatch - the maximum number of loci allowed to mistmatch between parent and offspring
#opt - has two options for output
#   1) parents and offspring with genotyeps, and number of mismatching loci
#   2) parents, offspring, and number of mismatching loci

exclusion <-function(adults, offspring, mismatch=0, opt= 1){

  # making sure MM is numeric and saving adults and offspring to other objects
  mismatch <- as.numeric(mismatch)
  Adults1 <- adults
  Offspring1 <- offspring

  #getting just the number of columns so gts and ids can be seperated
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]

  #getting the names for adults and offspring
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]

  #getting the number of genotypes, and how many adults and juveniles there are
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])

  #making all possible comparisons between adults and offspring using generic numbers
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  G

  #setting up two dfs with the gts of interest in the appropriate order
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
  head(Ads)
  head(Offs)

  #getting the IDs
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]
  IdnamesA
  IdnamesO

  #writing names to file
  write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)

  #reading them back in and combining them
  IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  Names

  #making a function to identify all possible matches between parent and offspring pairs (POPs) at a given locus
  matches <- function(matches){
    #comparing the genotypes of adult to offspring
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]

    #multipling them together and getting ride of NAs
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    #replacing those greater than 0 with 1, this is to ID those that do not match
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)

    #writing this to file and appending, I see the trick to avoid for loops
    write.table(f,file="Sort.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=T)
  }

  #getting the columns to apply this to
  z <- ncol(Ads)

  #applying the match function to all loci
  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)

  #reading results back in
  Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Observed

  #converting from long form to wide form
  a <- unique(Observed[,1])
  i<- NULL
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  } # end of wide form convert

  #adding up the riws for each column
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  head(Sorted)

  #identifying the matches
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)

  #identifying the locations of the names where matches occurred
  IDS <- which(stuff<(mismatch+1))

  #getting those ids for offspring and parents
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]

  #getting the number of putative pairs and filling a DF with them
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  Putativepairs

  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  PAdults
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  POffspring
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"

  #getting the gts for each parent and offspring drawing from f, which appeared earlier
  sorts <- function(sorts) {
    tell <- rbind(PAdults[f,],POffspring[f,])
    write.table(tell,file="Output_genotypes_tmp.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
  } #end of sort loop

  f <- length(PAdults[,1])

  #applying to all parents
  for(f in 1:f) lapply(f,sorts)

  #removing sort IdnamesA and IdnaemsO
  unlink("Sort.txt")
  unlink("IdnamesA.txt")
  unlink("IdnamesO.txt")

  #Outputfile processing==========================================================#
  putative<- read.table("Output_genotypes_tmp.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
  putative
  unlink("Output_genotypes_tmp.txt")
  if(opt == 1){

    return(putative)

  } #end option 1

  if(opt == 2){

    #seperating parents and offspring again
    parents <- putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
    offspring <- putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]
    parents
    offspring

    #now making two columns with IDs for parents and offspring and their number of MMs
    out1 <- cbind(parents,offspring)
    out2 <- out1[order(out1[,1]),-2]
    colnames(out2) <- c("Parent","Offspring","MM")
    head(out2)

    return(out2)
  } #end option 2

} #end of exclusion
