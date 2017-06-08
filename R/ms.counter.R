#################
### ms.counter
#################

#Description
#This fuction counts the mating success (MS) for each individual in that was used in parentage program to assign parents to offspring.

#Input Parameters:
#ped - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs
#sup.info - a dataframe with the first column as sample name and the second it's sex (Female ("F") or Male ("M")) and any number of columns after it

ms.counter <- function(sup.info, ped, missing.parent="UK"){

  mp <- missing.parent
  sup.info$mates <- 0
  sup.info$mates.uk <- 0

  for(i in 1:nrow(sup.info)){
    tmp1 <- NULL
    if(sup.info$sex[i] == "F"){

      tmp1 <- ped[ped[,2] == sup.info[i,1],]
      sup.info$mates[i] <- length(unique(tmp1[,3]))
      sup.info$mates.uk[i] <- length(tmp1[tmp1[,3] == mp, 3])
    }

    if(sup.info$sex[i] == "M"){
      tmp1 <- ped[ped[,3] == sup.info[i,1],]
      sup.info$mates[i] <- length(unique(tmp1[,2]))
      sup.info$mates.uk[i] <- length(tmp1[tmp1[,2] == mp, 2])
    }
  }
  return(sup.info)
} # end of ms.counter
