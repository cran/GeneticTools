# Version: 30-11-2012, Daniel Fischer

# This function uses the directional test/MW for testing for eQTL
eqtlDir <- function(genoGroups,gex){

  output <- c()
  for(i in 1:ncol(genoGroups))
  {
    if(sum(is.na(genoGroups[,i]))<(nrow(genoGroups))){
      PI <- length(table(genoGroups[,i]))
      # All individuals have the same genotype, do nothing
      if(PI==1){
	output[i] <- 1.25
      
      # Then the 2 groups comparison
      } else if (PI==2){
	output[i] <- gmw(gex[!is.na(genoGroups[,i])],genoGroups[!is.na(genoGroups[,i]),i],test="mw",type="permut",alternative="two.sided")

      # And the three group comparison
      } else if (PI==3){
	p1 <- gmw(gex[!is.na(genoGroups[,i])],genoGroups[!is.na(genoGroups[,i]),i],test="triple",type="permutation",alternative="greater",alg="Csub")
	p2 <- gmw(gex[!is.na(genoGroups[,i])],createGroups(genoGroups[!is.na(genoGroups[,i]),i],c(2,1,0)),test="triple",type="permutation",alternative="greater",alg="Csub")
	output[i] <- min(2*min(p1,p2),1)
      }
    } else {
      output[i] <- 1.5
    }
  }
  output
}