eqtlPlot <- function(snp=NULL, gene=NULL, eqtl=NULL,geno=NULL,gex=NULL, genoSamples=NULL){

  ## THIS DOESN*T WORK YET PROPERLY!!!
  if(!is.null(eqtl)){ 
    gex <- as.vector(t(eqtl$gex[,colnames(eqtl$gex)==gene]))
    names(gex) <- rownames(eqtl$gex)
    geno <- eqtl$geno
    genoSamples <- eqtl$genoSamples
  } else {
    gex <- gex
    geno <- geno
    genoSamples <- genoSamples
  }

  # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
  if(is.character(geno)==TRUE)
  {
    cat("Start reading the genotype information at",date(),"\n")
    genotData <- read.pedfile(file=paste(geno,".ped",sep=""),snps=paste(geno,".map",sep=""))
  } else {
    genotData <- geno
  }

  # Reducing the expression data to those rows, where we have also genotype information available
  gex <- gex[is.element(names(gex),genoSamples)]

  colOI <- which((geno$map$snp.names==snp)==TRUE)

  genoGroups <- as.numeric(as.vector(geno$genotypes[,colOI]))
  genoGroups <- rearrange(genoGroups,names(gex),genoSamples)

  values1 <- gex[genoGroups==1]
  values2 <- gex[genoGroups==2]
  values3 <- gex[genoGroups==3]

        PI <- length(table(genoGroups))
      # All individuals have the same genotype, do nothing
      if(PI==1){
	p <- 1.25
      
      # Then the 2 groups comparison
      } else if (PI==2){
	p <- gmw(gex[!is.na(genoGroups)],genoGroups[!is.na(genoGroups)],test="mw",type="permut",alternative="two.sided")

      # And the three group comparison
      } else if (PI==3){
	p1 <- gmw(gex[!is.na(genoGroups)],genoGroups[!is.na(genoGroups)],test="triple",type="permutation",alternative="greater",alg="Csub")
	p2 <- gmw(gex[!is.na(genoGroups)],createGroups(genoGroups[!is.na(genoGroups)],c(3,2,1)),test="triple",type="permutation",alternative="greater",alg="Csub")
	p <- min(2*min(p1,p2),1)
      }

  res <- list(v1=values1,v2=values2,v3=values3)
  boxplot(values1,xlim=c(0.5,3.5),at=1,ylim=c(min(gex),max(gex)),xaxt="n")
  boxplot(values2,at=2,add=TRUE)
  boxplot(values3,at=3,add=TRUE)
  axis(1,at=c(1,2,3),labels=c("A/A","A/B","B/B"))
  title(paste(snp,", p-value:",p))
}