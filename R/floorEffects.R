floorEffects_f1 <- function(x,y=0.95){
  y-ppois(1,x)
}


floorEffects <- function(confThresholds=c(0.95,0.99,0.999),popSize=121412){
  output <- matrix(nrow=length(confThresholds),ncol=2) %>%
  data.frame
  names(output) <- c("confidence","alleleFrequency")
  output$confidence <- confThresholds
  
  for(i in seq(along=output$confidence)){
    myMin <- uniroot(floorEffects_f1,c(0,1),output$confidence[i])$root
    output$alleleFrequency[i] <- myMin/popSize
  }
  return(output)  
}

# floorEffects() is looking for the AF at which a single allele is outside the specified poisson confidence, across a range of confidences