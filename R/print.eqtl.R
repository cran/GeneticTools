# Version: 30-11-2012, Daniel Fischer

`print.eqtl` <- function(x,which=NULL,sig=0.01,...){
  x <- x$eqtl
  xx <- x
  if(!is.null(which)==TRUE) xx <- x[which]

  X <- list()
  
  for(gene in 1:length(xx))
  {
      temp <- xx[[gene]]
      Nlocs <- length(table(temp$ProbeLoc))
      for(sub in 1:Nlocs)
      {
	temp <- xx[[gene]]$TestedSNP[xx[[gene]]$p.values<=sig,c(1,2,4,5,6)]
	temp2 <- xx[[gene]]$p.values[xx[[gene]]$p.values<=sig]
	temp <- cbind(temp,temp2)
	X[[gene]] <- temp
	colnames(X[[gene]]) <- c("Chr","SNP","Position","Allele1","Allele2","p.value")
      }
  }
  if(is.null(which)) which <- 1:length(x)
  names(X)<- names(x[which])
  print(X,...)
} 
