##########################
### randomization test ###
##########################

#Description
# This fuction was  written by Nick Sard in 2016. It is designed to take a two column
# data frame with reproductive success estimates (column 1) from two groups (column 2) and determine
# if the mean or variance are statistically significant based on a distrubution of n randomized comparisons.

randomization.test <- function(tmp,test="Two_Sided", paramater = "mean", group,n=1000){
  
  #getting the groups
  groups <- unique(tmp[,2])
  if(length(groups)!=2){stop("You need two groups to perform this function")}
  
  #defining groups
  g1 <- group
  g2 <- groups[groups!=group]
  
  #getting the values and seperating them by group
  vals <- tmp[,1]
  n.all <- length(vals)
  vals1 <- tmp[tmp[,2]==g1,1]
  vals2 <- tmp[tmp[,2]==g2,1]
  
  #getting counts of each group
  n1 <- length(vals1)
  n2 <- length(vals2)
  
  #doing the stats based on means
  if(paramater == "mean"){
    
    #calculating means for both groups and getting the difference
    mean1 <- mean(vals1)
    mean2 <- mean(vals2)
    real.diff <- mean1-mean2
    
    #saving this information for output df
    p1 <- mean1
    p2 <- mean2
    
    #doing the simulations
    sim.diffs <- NULL
    for(i in 1:n){
      sim.vals <- sample(x = vals,replace = F,size = n.all)
      sim.diff1 <- mean(sim.vals[1:n1])-mean(sim.vals[(n.all-n2+1):n.all])  
      sim.diffs <- c(sim.diffs,sim.diff1)
    } # end permutation
  } #end means if statement
  
  #doing the stats based on variance
  if(paramater == "variance"){
    
    #calculating vars for both groups and getting the difference
    var1 <- var(vals1)
    var2 <- var(vals2)
    real.diff <- var1-var2
    
    #saving this information for output df
    p1 <- var1
    p2 <- var2
    
    #doing the simulations
    sim.diffs <- NULL
    for(i in 1:n){
      sim.vals <- sample(x = vals,replace = F,size = n.all)
      sim.diff1 <- var(sim.vals[1:n1])-var(sim.vals[(n.all-n2+1):n.all])  
      sim.diffs <- c(sim.diffs,sim.diff1)
    } # end permutation
  } #end vars if statement
  
  #calculating p-value by comparing real.diff to simulated diffs
  if(test == "Two_Sided"){p.val <- length(sim.diffs[abs(sim.diffs)>abs(real.diff)])/n}
  if(test == "Less_Than"){p.val <- length(sim.diffs[sim.diffs<real.diff])/n}
  if(test == "Greater_Than"){p.val <- length(sim.diffs[sim.diffs>real.diff])/n}
  
  #putting all information together in a nice output format
  output <- data.frame(Group1 = g1, Group2 = g2,
                       Par1 = p1, Par2 = p2,
                       Diff = real.diff, Pval = p.val,stringsAsFactors = F)
  return(output)
} #end of permutation function