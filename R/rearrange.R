# Version: 21-12-2012, Daniel Fischer

# Changes:
# 21-12-2012
# genotypes can now be a vector or a matrix, needed for the eqtlPlot function

# This function takes a genotype matrix and rearrages its rows in that way, that the
# order is similar to a given expression data matrix.
# All rownames rows.gex should be included in genoSamples at this stage!!!
rearrange <- function(genoGroups,rows.gex,genoSamples){
  newOrder <- c()
  for(i in 1:length(rows.gex))
  {
    ## WARNING!!! IN CASE OF MULTIPLE MEASURES IN THE GENOTYPE, WE PICK ALWAYS THE FIRST OCCURENCE!!!!
    newOrder[i] <- which((rows.gex[i]==genoSamples)==TRUE)[1]
  }

  ifelse(is.matrix(genoGroups) , output <- genoGroups[newOrder,], output <- genoGroups[newOrder])
  ifelse(is.matrix(genoGroups) , rownames(output) <- genoSamples[newOrder], names(output) <- genoSamples[newOrder])
  
  output
}