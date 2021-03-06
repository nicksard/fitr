##################
### no.par.bayes()
##################

#Description
#This funciton takes a data.frame with parents and a data.frame with offspring and identify all parent offspring
#relationships with up to the number of mismatches bayes on bayesian priors and gives posterior probabilities for each
#take from SOLOMON - Mark Christie's library

#Input Parameters:
#adults - contains the sample names and genotypes of the parents
#offspring - contains the sample names and genotypes of the offspring
#reps - set the number of simulated datasets desired
#Ngts - set the number of simulated genotypes desired

no.par.bayes <- function(adults,offspring,reps=1000,Ngts=50000000)     {

  #loading in adults and offspring genotypes
  Adults1 <- adults
  Offspring1 <- offspring

  #making sure the number of reps is numeric and counting the number of loci
  reps <- as.numeric(reps)
  loci <- ncol(Adults1)

  #removing ID column
  Adults <- Adults1[,c(2:loci)]
  Offspring <- Offspring1[,c(2:loci)]

  #getting the number loci, again for Progress bar
  total <- ncol(Adults)/2

  #putting in a progress bar


  #Begin Master simulation function
  afreqs <- function(afreqs){

    #making locus names? *** need to figure out ** looks like used in progress bar
    locus_name <- L

    #setting up progress bar again


    #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
    vect <- c(Adults[,L],Adults[,L+1])
    alleles <- data.frame(table(vect))
    alleles <- alleles[order(alleles[,1]),]
    if (as.numeric(as.character(alleles[1,1]))==0) {alleles <- alleles[-1,]}           #removes missing data
    if(length(alleles[,1])==1) {                                                    #deals with monomorphic loci by adding 1 very strange allele (799)
      alleles <- cbind(vect[1],alleles[2])
      alleles[2,]<-c(799,1)}
    alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
    homos <- (alleles2[,3])^2                                                          #create homozygote allele frequencies
    homos2 <- cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
    hets <- t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
    hetfreq <- 2*(hets[,1]*hets[,2])
    hetvals <- t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
    hets2 <- cbind(hetvals,hetfreq)
    gfreqs <- rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
    csum <- cumsum(as.numeric(gfreqs[,3]))
    gfreqs1 <- cbind(gfreqs,csum)
    Nadults <- length(Adults[,1])
    Noffs <- length(Offspring[,1])


    #===============================================================================#end locus-specific HWE genotype frequency calculations
    alength <- length(alleles2[,1])
    for (Y in 1:reps) {
      positions <- 1:length(gfreqs[,1])
      sg3 <- sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                      #sample the repeated genotype positions, by the number of adults
      sadults <- gfreqs[sg3,1:2]                                                         #index gfreqs to create genotypes
      og3 <- sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                        #create juvenile genotyes
      soffs <- gfreqs[og3,1:2]
      soffs <- cbind(as.numeric(soffs[,1]),as.numeric(soffs[,2]))
      asp <- cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
      asp <- rbind(cbind(locus_name,0,0),asp)
      for (i in 1:Nadults) {
        parent1 <- as.numeric(sadults[i,1])                                                #first allele in parent
        parent2 <- as.numeric(sadults[i,2])                                                #second allele in parent
        p1 <- soffs-parent1
        p2 <- soffs-parent2
        pp1 <- which(p1[,1]==0)
        pp2 <- which(p1[,2]==0)
        allele1 <- unique(c(pp1,pp2))
        p21 <- which(p2[,1]==0)
        p22 <- which(p2[,2]==0)
        allele2 <- unique(c(p21,p22))
        Out51 <- cbind(parent1,length(allele1))
        Out52 <- cbind(parent2,length(allele2))
        Out53 <- cbind(0,Noffs-length(unique(c(allele1,allele2))))
        Out5 <- rbind(Out51,Out52,Out53)
        if(parent2==parent1) {Out5 <- Out5[-1,]}                                           #remove 1 of alleles count if homozygous
        if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
          diffs <- sum(Out5[,2])-Noffs                                                     #comment out to be more conservative!
          maxa <- max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
          pos <- which(Out5[,2]==maxa)
          Out5[pos,2]<-Out5[pos,2]-diffs}
        m1 <- match(Out5[,1],asp[,2])
        m2 <- asp[m1,3]+as.numeric(Out5[,2])
        asp[m1,3]<-m2
        asp<-asp
      }
      write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  L <- ncol(Adults)
  C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)

  #Bayes P-value==================================================================#
  OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  locname <- unique(OUT[,1])                                                         #compile calculations for each locus
  OUTALL <- NULL
  for (z in locname) {
    Loc1 <- OUT[which(OUT[,1]==z),]
    allfreqs <- unique(Loc1[,2])
    OUT2 <- NULL
    for (x in allfreqs) {
      a1<-Loc1[which(Loc1[,2]==x),]
      a2 <- sum(a1[,3])
      a3 <- cbind(x,a2)
      OUT2 <- rbind(OUT2, a3)
    }
    OUT3 <- cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
    OUTALL <- rbind(OUTALL, cbind(z,OUT3))
  }
  #Create multilocus genotypes, 50000 at a time, calculate number of loci that mismatch, and calculate freqs of shared alleles
  NL <- length(unique(OUTALL[,1]))
  ngtypes <- 1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
  Ngts <- as.numeric(Ngts)
  inreps <- 50000     #was tested as 10000 for SNPS
  repnumber <- round(Ngts/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  asp <- cbind(0:NL,rep(0,length(0:NL)))

  for (n in 1:repnumber) {

    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    distm <- apply(OUT, 1, function(x)sum(x == 1))
    distm2 <- data.frame(table(distm))
    m1 <- match(distm2[,1],asp[,1])
    m2 <- asp[m1,2]+distm2[,2]
    asp[m1,2]<-m2
    asp<-asp
  }

  #tabulate number of multilocus genotypes with x mismatches
  d2 <- asp
  d3 <- cbind(d2,d2[,2]/sum(d2[,2]))
  #===============================================================================# Create plot of exclusionary power.  Also, some necessary data formatting


  Adults<- adults
  Offs<- offspring
  Nads <- length(Adults[,1])
  Noffs <- length(Offs[,1])
  asp <- d3
  asp <- cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
  asp <- cbind(asp,cumsum(asp[,2]))
  asp <- cbind(asp,asp[,5]/max(asp[,5]))
  #find minimum Nloci to mismatch (could modify this) and perform exclusion with decided-upon mismatches==#
  distm <- cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
  #mismatch=min(which(round(distm[,6],1)==.9))                                    #deprecated
  mismatch <- which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
  Adults1 <- Adults                                                                   #begin exclusion
  Offspring1 <- Offs
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]
  write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  matches <- function(matches)
  {
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)
    write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  }
  z <- ncol(Ads)
  C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
  Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)
  IDS <- which(stuff<(mismatch+1))
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"
  sorts <- function(sorts)
  {
    tell <- rbind(PAdults[f,],POffspring[f,])
    write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  }
  f <- length(PAdults[,1])
  C1 <-  for(f in 1:f) lapply(f,sorts)
  unlink("Sort.txt")
  unlink("IdnamesA.txt")
  unlink("IdnamesO.txt")
  #Calculate phi for each number mismatching loci=================================#

  Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  observed <- data.frame(table(Putative[,2]))
  observed <- cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
  zerom <- 0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
  zerom2 <- which(is.na(match(zerom,observed[,1])))
  if (length(zerom2>0)) {observed <- observed[,3]
  for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}  #not really 0, to prevent divide by 0
  }   else {observed=observed[,3]}
  expected <- distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
  #expected=distm[1:(mismatch+1),3]*Nads*Noffs                                     #not using cumulative sum
  phi <- expected/observed
  phi <- replace(phi,which(phi>=1),1)
  Offs<- offspring
  actualTrue <- length(grep("Off",Offs[,1]))
  info <- cbind(actualTrue,expected,observed,phi)
  #calculate phi and index values ================================================#

  phibase <- phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
  observed <- observed[min(which(phi[]==1))-1]                                       #do the same for observed
  testob <- which(observed==0.000000001)
  phi2 <- cbind(1:length(phi),phi)
  if (length(testob>0)) {phi4 <- phi2[-testob,]} else {phi4 <- phi2}
  nmismatch <- min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
  index <- phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
  index <- index[which(index>-1)]
  if (length(index)>1) {
    if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
  phi <- phi[index]
  index <- index-1
  #Create Plot ===================================================================#
  pdf(file <- "Output_Dataset_Power.pdf")
  x <- 0:(length(info[,1])-1)
  y1 <- info[,3]
  y2 <- info[,2]
  p1 <- which(y2==0)
  y2 <- replace(y2,p1,.000000001)
  p1 <- which(y1<y2)
  y3 <- replace(y1,p1,y2[p1])
  y2 <- y2+1
  y3 <- y3+1
  par(mar = c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
  ats <- c(0,1,2,3,4,5,6,7,8,9)
  labs <- c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
  axis(side=2,ats,labs)
  lines(x,log10(y3),lwd=2)
  lines(x,log10(y2),lwd=2)
  points(x,log10(y3),pch=21,bg="green",cex=2)
  points(x,log10(y2),pch=21,bg="blue",cex=2)
  legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
  yphi <- y2/y3
  par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
  lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
  points(x,yphi,cex=2,pch=21,bg="gray")
  mtext("Number of Mismatching Loci",side=1,line=1.94)
  dev.off()
  info2 <- cbind(x,info[,4])
  colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
  write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  #===============================================================================#
  ngtypes <- 20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
  inreps <- 10000
  repnumber <- round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  #writes values at all loci, to be analyzed further below
  for (n in 1:repnumber) {
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      alleles3 <- alleles3[-which(alleles3[,2]==0),]
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    for (i in index) {                                                              #loop over numbers of mismatched loci
      if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                    #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distp <- cbind(i,DIST)
      write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }

  #Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
  OUT9 <- NULL
  for (n in index){
    Putative2 <- Putative[which(Putative[,2]==n),]
    write.table(Putative2, file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                                #need at least 1000 values , else assigned 0 (may be redundant now)
    OBS <- NULL
    afreqs <- function(afreqs) {
      PutL <- Putative2[,c(L+2,L+3)]
      PutLadults <- PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
      PutLoffs <- PutL[seq(from=2,to=length(PutL[,1]),by=2),]
      Puts <- cbind(PutLadults,PutLoffs)
      c1 <- Puts[,1]-Puts[,3]                                                    #find the matching alleles
      c2 <- Puts[,1]-Puts[,4]
      c3 <- Puts[,2]-Puts[,3]
      c4 <- Puts[,2]-Puts[,4]
      Puts2 <- cbind(Puts,c1,c2,c3,c4)
      P5 <- replace(Puts2[,5],which(Puts2[,5]!=0),-10)
      P6 <- replace(Puts2[,6],which(Puts2[,6]!=0),-10)
      P7 <- replace(Puts2[,7],which(Puts2[,7]!=0),-10)
      P8 <- replace(Puts2[,8],which(Puts2[,8]!=0),-10)
      P5 <- replace(P5,which(P5==0),1)
      P6 <- replace(P6,which(P6==0),1)
      P7 <- replace(P7,which(P7==0),1)
      P8 <- replace(P8,which(P8==0),1)
      Puts3 <- cbind(Puts,P5,P6,P7,P8)
      Puts4 <- cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
      alleles2 <- OUTALL[which(OUTALL[,1]==L),]
      alfreq1 <- alleles2[match(Puts4[,1],alleles2[,2]),4]
      alfreq2 <- alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
      alfreq3 <- alleles2[match(Puts4[,3],alleles2[,2]),4]
      alfreq4 <- alleles2[match(Puts4[,4],alleles2[,2]),4]
      Puts5 <- cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
      R1 <- replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
      R2 <- replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
      R3 <- replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
      R4 <- replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
      Puts6 <- cbind(R1,R2,R3,R4)
      Put_share <- apply(Puts6, 1, min)                                          #find row minimum
      Put_share2 <- apply(Puts4, 1, max)                                         #find shared allele name
      Put_share3 <- c(Put_share,Put_share2)
      OBS <<- cbind(OBS,Put_share3)
    }
    L <- ncol(Putative2)-2
    C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths <- length(OBS[,1])/2
    if (lengths==1) {OBA3 <- t(OBS[(lengths+1):(2*lengths),])} else {OBA3 <- OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file="True_shared_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS <- t(OBS[1:lengths,])} else {OBS <- OBS[1:lengths,]}   #formatting for if there is only a single pair
    obsp <- apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
    OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
  }
  #calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#

  OBA3<- read.table("True_shared_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
      vect <- c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
      alleles <- data.frame(table(vect))
      alleles <- alleles[order(alleles[,1]),]
      if (as.numeric(as.character(alleles[1,1]))<=0) {alleles <- alleles[-1,]}
      alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))
      gtrue <- sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
      OUT <- cbind(OUT,gtrue)

      if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
        mm <- alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
        write.table(cbind(r,mm),file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
      }
      }
    }
    for (i in index) {
      if (i==0) {DIST <- as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                   #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distt<-cbind(i,DIST)
      write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  #Calculate lamdaphi=============================================================#

  Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Putadults <- Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
  Putoffs <- Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
  Names <- cbind(as.character(Putadults),as.character(Putoffs))
  empirical <- cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
  distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  P2 <- NULL
  for (i in index) {                                                              #loop over numbers of mismatched loci
    empirical2 <- empirical[which(as.numeric(as.character(empirical[,3]))==i),]
    if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
    if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
      a3 <- distp[which(distp[,1]==i),2]
      DIST<-as.data.frame(a3)
      P <- NULL
      for (b in 1:length(empirical2[,1])) {
        p1 <- length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
        if (p1==0) {p1=0.00001}
        p2 <- cbind(empirical2[b,1],empirical2[b,2],p1)
        p3 <- cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
        P <- rbind(P,p3)
      }
    }
    P2<-rbind(P2,cbind(i,P))
  }
  lamdaphi <- as.numeric(as.character(P2[,5]))
  #Calculate lamda|phic ==========================================================#

  lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  mm1 <- which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
  Observed[mm1,2] <- 1                                                            #replace NAs with 1
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  lamdaphic <- apply(U, 1, prod)
  l1 <- length(which(OUT9[,2]==0))
  if (l1>0) lamdaphic <- c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
  P3 <- cbind(P2,lamdaphic)
  for (i in index) {                                                              #loop over numbers of mismatched loci
    e2 <- P3[which(as.numeric(as.character(P3[,1]))==i),]
    if(length(e2)==6) {e2 <- t(e2)}                                                 #deals with one indiviudal formatting
    if (e2[1,5]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
      a3 <- lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
      DIST<-as.data.frame(a3)
      for (b in 1:length(e2[,1])) {                                                 #calculate p values
        p1 <- length(which(DIST[,1] <= e2[b,6]))
        p2 <- p1/length(DIST[,1])
        write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      }
    }
  }
  lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  #Put it all together with Bayes theorem!========================================#

  vals <- cbind(P2,lamdaphic[,1])
  philength <- cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
  phis <- rep(philength[,2],philength[,3])
  vals <- cbind(vals,phis)
  colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
  phi <- as.numeric(as.character(vals[,7]))
  lamdaphi <- as.numeric(as.character(vals[,5]))
  lamdaphic <- as.numeric(as.character(vals[,6]))
  lamdaphi <- replace(lamdaphi,which(lamdaphi==0),1)
  lamdaphic <- replace(lamdaphic,which(lamdaphic==0),1)
  pval <- (lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
  pval <- cbind(vals[,2],vals[,3],vals[,1],pval)
  pval <- pval[order(as.numeric(pval[,4])),]
  colnames(pval) <- c("Adult","Offspring","NL_mismatch","Probability of pair being false given frequencies of shared alleles")
  #write.table(vals, file="Posterior_Components_Bayes.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)

  unlink("False_allele_freqs.txt")                                                #clear all sims files
  unlink("True_allele_freqs.txt")
  unlink("Shared_allele_freqs.txt")
  unlink("lamdaphic.txt")
  unlink("Putative.txt")
  unlink("True_shared_freqs.txt")
  unlink("Output_genotypes.txt")
  unlink("*.sims")
  rm(list=ls())
} #end of exclusion
