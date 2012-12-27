# Version: 22-12-2012, Daniel Fischer

# Changes: 22-12-2012
# added the function plotIt() in order to avoid redundancies

`plot.eqtl` <- function(x, file=NULL, which=NULL, sig=0.01, ...){
  windowSize <- x$windowSize
  x <- x$eqtl
  
  if(!is.null(which)==TRUE) x <- x[which]

  if(is.null(file))
  {
    plotIt(x=x,sig=sig,windowSize=windowSize)
  } else {
    pdf(file=file,width=10,height=10)
      plotIt(x=x,sig=sig,windowSize=windowSize)
    dev.off()
  }
  invisible()
} 


plotIt <- function(x,sig,windowSize){
  for(gene in 1:length(x))
  {
      temp <- x[[gene]]
      Nlocs <- length(table(temp$ProbeLoc))
      for(sub in 1:Nlocs)
      {
	subPos <- temp$ProbeLoc==sub
	minX <- min(temp$TestedSNP[subPos,4])
	maxX <- max(temp$TestedSNP[subPos,4])
	plot(c(-10,-10),ylim=c(0,1.5),xlim=c(0,20),xlab="Chromosomal Position in MB",ylab="p-value",main=paste(names(x)[gene],"-",sub),sub=paste("Chr",temp$GeneInfo[sub,1],":",minX,"-",maxX),yaxt="n",xaxt="n")
	axis(2,at=c(seq(0,1,0.2),1.25,1.5),labels=c(seq(0,1,0.2),"Homoc.","NA"))
	axis(1,at=seq(0,20,4),labels=seq(round(minX/10^6,1),round(maxX/10^6,1),length.out=6))

	xPos <- (temp$TestedSNP[subPos,4]-minX)/(windowSize*10^5)
	yPos <- temp$p.values[subPos]
	# eliminate rounding errors (drop all p-values which are a bit larger than 1 to 1)
	tempPos <- (yPos>1) & (yPos < 1.2)
	yPos[tempPos] <- 1
	col <- rep("green",length(subPos))
	col[temp$p.values[subPos] <= sig] <- "red" 
	col[yPos[subPos] > 1] <- "black"

	points(xPos,yPos,col=col,pch=20)
	
	lines(c(-10,20),c(1,1),lty="dotted")
	lines(c(-10,20),c(sig,sig),lty="dotted")

	# Plot the gene position
	xG <- (mean(unlist(c(temp$GeneInfo[sub,][2:3])))-minX)/(windowSize*10^5)
	lines(c(xG,xG),c(-1,2),lty="dashed")
      }
  }
}