# Version: 22-12-2012, Daniel Fischer

eQTL <- function(gex, geno, xAnnot, xSamples=NULL, genoSamples, windowSize=0.5, method="LM"){
  
  # If the annotations are given as data frame, we will transform them into a list
  if(is.data.frame(xAnnot)){
    xAnnot <- makeAnnotList(xAnnot)
  }

  # Input checks
  single <- FALSE
  th <- windowSize
  method <- match.arg(method,c("LM","directional"))

  # If only one gene is given it could be a vector, transform it here to a column matrix
  if(is.vector(gex)){
    single <- TRUE
    tempNames <- names(gex)
    if(!is.null(xSamples)) tempNames <- xSamples
    gex <- matrix(gex,ncol=1)
    rownames(gex) <- tempNames
    colnames(gex) <- names(xAnnot)
    gexColNames <- colnames(gex)
  }

  # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
  if(is.character(geno)==TRUE)
  {
    cat("Start reading the genotype information at",date(),"\n")
    genotData <- read.pedfile(file=paste(geno,".ped",sep=""),snps=paste(geno,".map",sep=""))
  } else {
    genotData <- geno
  }

  # Sample statistics
  overlap <- is.element(rownames(gex),genoSamples)
  olPerc <- round(sum(overlap)/nrow(gex)*100,1)
  if(olPerc==0) stop("No matching expression / genotype sample names!\n")
  cat("We have for",olPerc,"% of the samples in the expression data the genotype information. \n")

  # Location statistics
  overlap <- is.element(colnames(gex),names(xAnnot))
  olPerc <- round(sum(overlap)/ncol(gex)*100,1)
  if(olPerc==0) stop("No matching expression probe names / probe name annotations!\n")
  cat("We have for",olPerc,"% of the expression data probes the probe annotations. \n")

  # Probe statistics
  matchingProbes <- colnames(gex)[is.element(colnames(gex),names(xAnnot))]
  cat("We will investigate for",length(matchingProbes),"probes possible eQTLs! \n")
  result <- list()

  # Reducing the expression data to those rows, where we have also genotype information available
  gex <- gex[is.element(rownames(gex),genoSamples),]
  if(single==TRUE){
    gex <- t(t(gex))
    colnames(gex) <- gexColNames
  }

  # Now go through all possible probes
  eqtl <- list()
  for(probeRun in 1:length(matchingProbes))
  {
    # Do that for each possible location of the probe (might not be unique...)
    tempAnnot <- xAnnot[[which((names(xAnnot)==matchingProbes[probeRun])==TRUE)]]
    eqtlTemp <- list()
    for(tempRun in 1:nrow(tempAnnot))
    {
      # Temporary values
      SNPloc <- getSNPlocations(genotData$map,tempAnnot[tempRun,],th=th)
      SNPmatrix <- genotData$genotypes[,SNPloc$SNPcol]
      genoGroups <- as(SNPmatrix,"numeric")
      genoGroups <- rearrange(genoGroups,rownames(gex),genoSamples)
      
      # Type of eQTL
      if(method=="LM"){
	eqtlTemp[[tempRun]] <- list(ProbeLoc=rep(tempRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlLM(genoGroups,gex[,probeRun]))
      } else if(method=="directional"){
	eqtlTemp[[tempRun]] <- list(ProbeLoc=rep(tempRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlDir(genoGroups,gex[,probeRun]))
      }
    }
    
    # Join the output
    eqtl[[probeRun]] <- joinEQTL(eqtlTemp)
    eqtl[[probeRun]]$GeneInfo <- tempAnnot
    if((probeRun %% 100) == 0) cat ("We managed to calculate",probeRun,"eQTLs at",date(),"\n")
  }

  # Return the result
  names(eqtl) <- matchingProbes
  result <- list(eqtl=eqtl,gex=gex, geno=geno, xAnnot=xAnnot, genoSamples=genoSamples, windowSize=windowSize)
  class(result) <- "eqtl"
  result
}



