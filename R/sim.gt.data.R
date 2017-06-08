#####################
### sim.gt.data() ###
#####################

#Description
#This fuction is a modified version of Mark Christie's script to create parents and offspring.
#Note that one can create random genotyping error, add unrelated indviduals, and add full siblings in the mix

#Input Parameters:
#afreq - this is a data frame with two rows in it, the first is the locus, and the second is the allele frequencies
#Nparents - defining thenumber of parents one wants
#Noffs_perpair - defining the number of offspring per parent one wants
#error - defining the amount of random error in the dataset
#Nunrelated - defining the number of unrelated individuals wanted in the dataset
#Nsibs - defining the number of full siblings one wants

#Create Simulated Data Sets==========================================#

sim.gt.data <-  function(afreq, Nparents= 1, Noffs_perpair= 1, error= 0, Nunrelated= 0, Nsibs=0, write2file = F){

  #making sure everything is numeric
  Nparents <- as.numeric(Nparents)
  Noffs_perpair <- as.numeric(Noffs_perpair)
  error <- as.numeric(error)
  Nunrelated <- as.numeric(Nunrelated)

  #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
  Nadults <- Nparents*2

  #getting afreqs input
  afreqs <- afreq
  afreqs

  #making a function to simulate genotypes for parents

  OUT <- NULL
  sims <- function(sims){

    #table allele frequencies
    #getting just the frequencies
    alleles2 <- afreqs[afreqs[,1] == loci[i],]
    alleles2

    #changing generic name to three character genotypes
    alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])
    alleles3

    #create homozygote allele frequencies
    homos <- (alleles3[,2])^2
    homos
    homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
    homos2

    #create heterozygote allele frequencies

    #first create all possible combinations of freqs
    hets <- t(combn(alleles3[,2],2))
    hets

    #now make expected freqs
    hetfreq <- 2*(hets[,1]*hets[,2])
    hetfreq

    #create heterozygote allele names
    #now getting the het. genotype combinations
    hetvals <- t(combn(as.character(alleles3[,1]),2))
    hetvals

    #combing them
    hets2 <- cbind(hetvals,hetfreq)
    hets2

    #combine hets and homos and create genotypes
    gfreqs <- rbind(hets2,homos2)
    gfreqs

    #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
    n <- 1000000

    #create genotypes(by coloumn, 1 for each allele)
    #for the first allele

    gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))
    gfreqs1

    #now the second
    gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
    gfreqs2

    #combining them
    gtypes <- cbind(gfreqs1,gfreqs2)
    head(gtypes)

    #mixing them up
    gtypes <- gtypes[sample(1:length(gtypes[,1]),replace= F),]
    head(gtypes)

    #now getting a random sample of the gentoypes for the parents
    sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
    sg1

    OUT<<-cbind(OUT,sg1)

  } # end of sim function

  #defining the number of loci there are
  loci.length <- length(unique(afreqs[,1]))
  loci <- unique(afreqs[,1])
  loci

  #doing the simulation process for each locus
  for(i in 1:length(loci)){
    lapply(1,sims)
  } # end for loop
  OUT
  warnings()
  OUT
  #saving OUT into object called parents
  parents <- OUT
  head(parents)

  #next making a matrix to add to parents so that gts in every two columns don't overlap
  #first making a matrix the same size of parents and filling it ever two columns with varying
  #numbers in the 1000s
  c <- c(1:(ncol(OUT)))
  odd <- 2*(unique(round(((c-2))/2)))+1
  l <- length(odd) * 1000
  codes <- seq(from= 1,to= l,by= 1000)
  cols <- sort(rep(codes,2))-1
  Anumbs <- matrix(data= cols,nrow= Nadults,ncol= ncol(OUT),byrow= T)
  Anumbs

  #now adding it to parents
  parents <- as.numeric(parents)+Anumbs
  parents

  #create full sib families (go down the list in pair) ======================#
  OUT2 <- NULL
  sims <- function(sims){

    #first grabbing the genotypes of the parents
    p1 <- parents[i,]
    p1
    p2 <- parents[i+1,]
    p2
    #defining the order of teh alleles
    als <- rep(1:2,length(p1)/2)
    als

    #defining the number of offspring that need to be full sibs
    Noffs <- Noffs_perpair
    Noffs

    #using a for loop to making full sib genotypes
    OUT2 <- NULL
    for (b in 1:Noffs){

      #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
      #sampling alleles from parent one
      pos1 <- sample(als,length(p1)/2,replace <- TRUE)
      pos1

      #sampling alleles from parent two
      pos2 <- sample(als,length(p1)/2,replace <- TRUE)
      pos2

      #getting the position for those alleles from parent one
      pos11 <- pos1+(seq(0,(length(p1)-1),2))
      pos11

      #getting the position for those alleles from parent two
      pos22 <- pos2+(seq(0,(length(p2)-1),2))
      pos22

      #getting those alleles from parent one
      o1 <- p1[pos11]
      o1
      #getting those alleles from parent two
      o2 <- p2[pos22]
      o2

      #putting them together to make the offsprings genotype
      o3 <- sort(c(o1,o2))
      o3

      #identifying what parents the offspring came fromm
      o3 <- t(c(i,o3))
      o3

      #writting to file and appending
      write.table(o3,file <- "SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
    } #end of full sib for loop
  } # end of new sims function

  #defining which parents are going to have full sibs
  my.picks <- seq(from=1,to=Nadults-1, by=2)
  my.picks

  #making full sib gts for each parent pair
  for(i in my.picks){
    lapply(i,sims)
  } # end of full sib genotyping loop

  #removing the 1000s from the gts
  Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
  parents <- as.numeric(parents)-Anumbs
  parents

  #getting parent genotypes
  #first the moms
  Moms <- parents[seq(from=2,to= Nadults,by= 2),]
  Moms
  #now the dads
  Dads <- parents[seq(from=1,to= Nadults,by= 2),]
  Dads

  #reading SimOffs.txt back in and removing the id column
  Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Offs  <-  as.matrix(Offs2[,-1])

  #removing issue with 1000s
  Anumbs <- matrix(cols,length(Offs[,1]),ncol(OUT),byrow = T)
  Offs <- as.numeric(Offs)-Anumbs
  Offs

  #now making offspring names
  Noffs_perpair
  if (Noffs_perpair>1){

    Offnames <- ceiling(Offs2[,1]/2)
    Offnames2 <- paste0(Offnames,".",1:length(Offnames))
    Offnames2
    Offs <- cbind(paste("Offspring",Offnames2),Offs)
    Offs

  } else {

    Offnames <- ceiling(Offs2[,1]/2)
    Offs <- cbind(paste("Offspring",Offnames),Offs)
  } # end of if statement about offspring

  Moms <- cbind(paste("Mom",1:length(Moms[,1])),Moms)
  Moms
  Dads <- cbind(paste("Dad",1:length(Dads[,1])),Dads)
  Dads

  #add error======================================================================#
  error
  #first for the dads

  #getting the gts
  Dadsg <- Dads[,-1]

  #figuring out how many gts to make errors in
  ldad <- length(Dadsg)*error
  ldad
  #randomly selecting where to put the errors
  pdad1 <- sample(1:length(Dadsg),ldad,replace= FALSE)
  pdad1

  #randomly selecting some gets to replace in the those locations
  pdad2 <- Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
  pdad2

  #actually putting in the error
  Dadsg2 <- replace(Dadsg,pdad1,pdad2)
  Dadsg2
  #replacing the old Dads gts
  Dads <- cbind(Dads[,1],Dadsg2)
  Dads

  #doing the same thing for the kids
  Offsg=Offs[,-1]
  loff=length(Offsg)*error
  poff1=sample(1:length(Offsg),loff,replace=FALSE)
  poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
  Offsg2=replace(Offsg,poff1,poff2)
  Offs=cbind(Offs[,1],Offsg2)

  #not sure why he does not add error in to the moms I guess that would double the error rate?
  #Momsg <- Moms[,-1]
  #ldad <- length(Momsg)*error
  #pdad1 <- sample(1:length(Momsg),ldad,replace=FALSE)
  #pdad2 <- Momsg[sample(1:length(Momsg),ldad,replace=FALSE)]
  #Momsg2 <- replace(Momsg,pdad1,pdad2)
  #Moms <- cbind(Moms[,1],Momsg2)

  #===============================================================================#
  #making generic locus names
  colsnmz <- rep("Locus",ncol(Offs)-1)
  colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
  colsnmz3 <- paste(colsnmz,colsnmz2)
  colsnmz3 <- c("IDs",colsnmz3)
  colsnmz3

  #usig them to make the colnames for the Offs, Dads, and Moms
  colnames(Offs)<- colsnmz3
  colnames(Dads)<- colsnmz3
  colnames(Moms)<- colsnmz3

  #Add unrealated individuals <- ====================================================#
  Nunrelated

  if (Nunrelated>0){

    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nunrelated*3
    afreqs <- afreq
    OUT <- NULL

    #making the allele frequencies again, this way they are completely random and unrelated
    #doing the same thing as lines 34-113
    sims <- function(sims) {
      alleles2 <- afreqs[which(afreqs[,1] == z),]
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])

      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)

      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))
      hetfreq <- 2*(hets[,1]*hets[,2])

      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))
      hets2 <- cbind(hetvals,hetfreq)

      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)

      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000

      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]

      OUT <<-cbind(OUT,sg1)

    } #end of sims function

    #now applying function
    z1 <- length(unique(afreqs[,1]))
    for(z in 1:z1) {lapply(z,sims)}

    #repeating lines 122-204
    parents <- OUT
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    parents <- as.numeric(parents)-Anumbs

    #to called this thing unrelated
    unrelated <- cbind(paste("Individual",1:length(parents[,1])),parents)
    Lid <- length(unrelated[,1])/3

    #splitting the unrelated individuals among dads, moms, and kids
    m1 <- unrelated[1:Lid,]
    d1 <- unrelated[(Lid+1):(Lid*2),]
    j1 <- unrelated[((Lid*2)+1):(Lid*3),]

    #now for the column names again
    colsnmz <- rep("Locus",ncol(Offs)-1)
    colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
    colsnmz3 <- paste(colsnmz,colsnmz2)
    colsnmz3 <- c("IDs",colsnmz3)
    colnames(Offs)<- colsnmz3
    colnames(Dads)<- colsnmz3
    colnames(Moms)<- colsnmz3

    #adding the unrelated individuals in the three files
    Moms <- rbind(Moms,m1)
    Dads <- rbind(Dads,d1)
    Offs <- rbind(Offs,j1)

  } #end of unrelated if statement

  #removing the file SimOff.txt from working directory
  unlink("SimOffs.txt")

  #Begin add siblings=============================================================#
  #This scripts creates a desired number of siblings and splits them between adults and offspring files
  Nsibs <- as.numeric(Nsibs)
  if (Nsibs > 0) {
    #could change this, but as stands only evaluating 1 pair of siblings per parent-pair
    Noffs_per_pair <- 2

    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nsibs*2
    afreqs <- afreq

    #repeating allele frequency simulation from above
    OUT <- NULL
    sims <- function(sims){
      alleles2 <- afreqs[which(afreqs[,1]==z),]

      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])

      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)

      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))
      hetfreq <- 2*(hets[,1]*hets[,2])

      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))
      hets2 <- cbind(hetvals,hetfreq)

      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)

      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000

      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]

      OUT <<-cbind(OUT,sg1)

    } #end of allele frequency simulation

    #applyingn to making genotypes of parengs
    z <- length(unique(afreqs[,1]))
    for(z in 1:z) {lapply(z,sims)}
    parents <- OUT

    #repeating the part where he addes 1000 to each line... still not sure why he does this
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs


    #create full sib families (go down the list in pair)============================#
    OUT2 <- NULL
    sims <- function(sims)  {
      N <- 1:length(parents[,1])
      u <- sample(N,1)
      u2 <- sample(N,1)
      p1 <- parents[u,]
      p2 <- parents[u2,]
      als <- rep(1:2,length(p1)/2)

      #number of offspring per pair
      Noffs <- Noffs_per_pair
      OUT2 <- NULL
      for (b in 1:Noffs){

        #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
        pos1 <- sample(als,length(p1)/2,replace=TRUE)
        pos2 <- sample(als,length(p1)/2,replace=TRUE)
        pos11 <- pos1+(seq(0,(length(p1)-1),2))
        pos22 <- pos2+(seq(0,(length(p2)-1),2))
        o1 <- p1[pos11]
        o2 <- p2[pos22]
        o3 <- sort(c(o1,o2))
        o3 <- t(c(z,o3))
        write.table(o3,file="SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
      } #end of for loop
    } # end of simulation for sibs

    #applying it to the parents
    z <- length(parents[,1])

    #used to move down list, now sample randomly so is just used to produce wanted number of offspring
    C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,sims)
    Dads <- parents[seq(from=1,to=length(parents[,1]),by=2),]
    Moms <- parents[seq(from=2,to=length(parents[,1]),by=2),]

    #see code before functions for adding 1000s (here am removing 1000s)
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)-Anumbs
    Dads <- parents[seq(from=1,to=length(parents[,1]),by=2),]
    Moms <- parents[seq(from=2,to=length(parents[,1]),by=2),]

    #reading in the kids again
    Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    Offs  <-  as.matrix(Offs2[,-1])
    Anumbs <- matrix(cols,length(Offs[,1]),ncol(OUT),byrow <- T)
    Offs <- as.numeric(Offs)-Anumbs

    #naming offspring
    if (Noffs_per_pair>1){
      #naming of offpsring
      Offnames <- ceiling(Offs2[,1]/2)
      Offnames2 <- paste(Offnames,".",1:length(Offnames))
      Offs <- cbind(paste("Sibling",Offnames2),Offs)
    } else {
      Offnames <- ceiling(Offs2[,1]/2)
      Offs <- cbind(paste("Sibling",Offnames),Offs)
    } # end of offspring naming if statement

    Moms <- cbind(paste("Mom",1:length(Moms[,1])),Moms)
    Dads <- cbind(paste("Dad",1:length(Dads[,1])),Dads)


    #add error======================================================================#
    Offsg <- Offs[,-1]
    loff <- length(Offsg)*error
    poff1 <- sample(1:length(Offsg),loff,replace=FALSE)
    poff2 <- Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
    Offsg2 <- replace(Offsg,poff1,poff2)
    Offs <- cbind(Offs[,1],Offsg2)

    #===============================================================================#
    colsnmz <- rep("Locus",ncol(Offs)-1)
    colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
    colsnmz3 <- paste(colsnmz,colsnmz2)
    colsnmz3 <- c("IDs",colsnmz3)
    colnames(Offs)<- colsnmz3

    #calculate shared alleles among pairs of siblings===============================#
    sib1 <- Offs[seq(from=1,to=length(Offs[,1]),by=2),]
    sib2 <- Offs[seq(from=2,to=length(Offs[,1]),by=2),]
    write.table(sib1,file="Sib1.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(sib2,file="Sib2.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)

    #appending siblings to file
    #getting dad siblings
    dadsibs <- read.table("Dads_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    sib1 <- read.table("Sib1.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    dadsibs <- rbind(dadsibs,sib1)
    #write.table(dadsibs,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)

    #putting them in the kid file
    juvsibs <- read.table("Juveniles_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    sib2 <- read.table("Sib2.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    juvsibs <- rbind(juvsibs,sib2)
    #write.table(juvsibs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)

    #adding them to the end of the parent files
    Dads <- rbind(Dads,dadsibs)
    Offs <- rbind(Offs,juvsibs)

    #removing files
    unlink("Sib1.txt")
    unlink("Sib2.txt")
    unlink("SimOffs.txt")
  } # end of sibling function


  if(write2file == T){
    print("Mom, Dad, and Offspring genotypes have been saved to current working directory")
    write.table(Dads,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(Moms,file="Moms_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(Offs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
  } else {

    Parents <- rbind(Moms,Dads)
    Parents <- data.frame("Parents",Parents,stringsAsFactors=F)
    colnames(Parents) <- c("Type", colnames(Parents)[-1])

    Offspring <- data.frame("Offspring",Offs,stringsAsFactors=F)
    colnames(Offspring) <- c("Type",colnames(Offspring)[-1])

    output <- data.frame(rbind(Parents,Offspring),stringsAsFactors=F)
    return(output)
  }
} # end of sim.gt.data function
