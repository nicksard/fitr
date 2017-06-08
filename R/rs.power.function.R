#####################
### rs.rand.power ###
#####################

#Description
# This fuction was  written by Nick Sard in 2016. It is designed to determine the power of 
# a randomization test with some given number of randomized replicates (reps)
# to detected a statistical difference between two samples (n1,n2) that come from
# two different negative binomial distrbutions that have mean1 and mean2 and variances (var1, var2).
#
# Power is defined as the probability of rejecting a false null hypothesis. 
# This function uses defined values of n1, n2, mean1, mean2, var1, var2 and
# conducts a randomization test n times, and calculates power as the 
# porportion of n randomization tests that detected a statistically significant
# difference between the two samples at some alpha level. This function is written
# with some flexibility so that it can compare means or variances, and two sided
# or one sided (great than or less than) alternative hypotheses.

rs.rand.power <- function(n1=30,n2=30,mean1=20,mean2=21,var1=21,var2=20,n=1000,reps=10000,stat="mean",type="Two_sided",alpha=0.05){
  
  #need to load library MASS because we are using a Negative Binomial distribution
  require(MASS)
  i <- NULL
  OUT <- NULL
  for(i in 1:n){
    #print(i)
    #simulating neg.bin distributed rs data for dataset 1
    if(((var1-mean1))==0){
      print("Mean and Variance are equal in sample 1. Setting theta to 0.001")
      theta=0.001
    } else {
      theta1 <- abs(mean1^2/(var1-mean1))
    } #end of sample1 if statement
  
    #simulating neg.bin distributed rs data for dataset 2
    if(((var2-mean2))==0){
      print("Mean and Variance are equal in sample 1. Setting theta to 0.001")
      theta=0.001
    } else {
      theta2 <- abs(mean2^2/(var2-mean2))
    } #end of sample2 if statement
    
    #simulatining two different rs distributions
    rs1 <- rnegbin(n = n1,mu = mean1,theta = theta1)
    rs2 <- rnegbin(n = n2,mu = mean2,theta = theta2)
    
    #combing for randomization
    rss <- c(rs1,rs2)
    
    if(stat == "mean"){
      
      #calculating the difference between them
      true.diff <- mean(rs1)-mean(rs2)
      
      #making a randomization function
      randomization <- function(){
        #shuffling gata and putting into two groups
        gs <- sample(x = rss,size = n1+n2,replace = F)
        g1 <- gs[1:n1]
        g2 <- gs[(n1+1):(n1+n2)]
        m1 <- mean(g1)
        m2 <- mean(g2)
        m.diff <- m1-m2
        return(m.diff)
      } #end of randomization function  
    } #end of mean if statement
    
    if(stat == "variance"){
      
      #calculating the difference between them
      true.diff <- var(rs1)-var(rs2)
      
      #making a randomization function
      randomization <- function(){
        #shuffling gata and putting into two groups
        gs <- sample(x = rss,size = n1+n2,replace = F)
        g1 <- gs[1:n1]
        g2 <- gs[(n1+1):(n1+n2)]
        v1 <- var(g1)
        v2 <- var(g2)
        v.diff <- v1-v2
        return(v.diff)
      } #end of randomization function  
    } #end of mean if statement
    
    reps1 <- sort(replicate(n = reps,randomization()))
    table(round(reps1))
    if(type == "Two_sided"){
      a.reps<- abs(reps1)
      a.true.diff <- abs(true.diff)
      pval <- length(a.reps[a.reps>a.true.diff])/n
    } #end of two sided if statement
    
    if(type == "One_sided"){
      pval <- length(reps1[reps1<true.diff])/n
    } #end of one sided if statement
    OUT<-c(OUT,pval)
  } #end of for loop

  #calculating power
  power <- length(OUT[OUT<=alpha])/n
  
  #creating a table to summarize the simulation and provide power estimate
  Info <- c("N1","N2","mean1","mean2","var1","var2","stat","type","alpha","Power")
  Values <- c(n1,n2,mean1,mean2,var1,var2,stat,type,alpha,power)
  output <- data.frame(Info,Values,stringsAsFactors = F)
  return(output)
} # end of power function

