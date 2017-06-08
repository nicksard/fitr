#loading my own functions
bootstrap.cis <- function(boot.values,point.est,ci=95){
  
  #first determine which values in my boot.values will be my picks for the ci
  orig.ci <- ci
  ci <- ((100-ci)/2)/100
  n.boot.values <- length(boot.values)
  ci.1 <- ci*n.boot.values
  ci.2 <- n.boot.values-ci.1
  
  #now calculating the variation around my point estimate and sorting from smallest to largest
  boot.values <- sort(boot.values-point.est)
  
  #now acutally getting my cis
  ci.1 <- boot.values[ci.1]
  ci.2 <- boot.values[ci.2]
  ci.low <- point.est-ci.2
  ci.high <- point.est-ci.1
  
  #now making a df to return
  df <- data.frame(point.est,ci.low,ci.high,ci.1,ci.2)
  colnames(df) <- c("Point.Estimate",paste0("CI.",orig.ci,".Low"),paste0("CI.",orig.ci,".High"),"ci1","ci2")
  df
}