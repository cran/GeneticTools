# Version: 30-11-2012, Daniel Fischer

getSNPlocations <- function(genotInfo,annot,th){
  th <- th * 10^6
  chrSNPs <- genotInfo[genotInfo[,1]==as.character(annot[1,1]),]
  lowSNPs <- chrSNPs[chrSNPs[,4]>(annot$Start-th),]
  SNPs <- lowSNPs[lowSNPs[,4]<(annot$End+th),]
  output <- list(SNPloc=SNPs,SNPcol=as.numeric(rownames(SNPs)))
  output
}