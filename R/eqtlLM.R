# Version 30-11-2012, Daniel Fischer

# This function applies a linear Regression model to the data and tests of the slope is equal to zero or not.
eqtlLM <- function(genoGroups,gex){

  output <- c()
  for(i in 1:ncol(genoGroups))
  {
    if(sum(is.na(genoGroups[,i]))<(nrow(genoGroups))){
      fitData <- data.frame(x=genoGroups[,i],y=gex)
      temp <- lm(y~x,data=fitData)
      ifelse(is.na(temp$coefficients[2]), output[i] <- 1.25 ,  output[i] <- summary(temp)$coefficients[2,4])
    } else {
      output[i] <- 1.5
    }
  }
  output
}