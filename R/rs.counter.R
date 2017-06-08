#################
### rs.counter
#################

#Description
#This fuction counts the reproductive success (RS) for each individual in that was used in parentage program to assign parents to offspring.

#Input Parameters:
#ped - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs
#sup.info - a dataframe with the first column as sample name and the second it's sex (Female ("F") or Male ("M")) and any number of columns after it

rs.counter <- function(sup.info, ped){

  #making a column for RS
  sup.info$rs <- 0

  #making sure the necessary columns are in character form
  ped[,2] <-as.character(ped[,2])
  ped[,3] <-as.character(ped[,3])
  sup.info[,1] <- as.character(sup.info[,1])

  j <- NULL
  #here is the for loop
  for(j in 1:nrow(sup.info)){
    if(sup.info[j,2] == "F"){
      sup.info$rs[j] <- nrow(ped[ped[,2] == sup.info[j,1],])
    }

    if(sup.info[j,2] == "M"){
      sup.info$rs[j] <- nrow(ped[ped[,3] == sup.info[j,1],])
    }
  }
  return(sup.info)
} #end of rs.counter
