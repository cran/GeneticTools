# Version: 30-11-2012, Daniel Fischer

summary.eqtl <- function(object, ...){
  
  cat("EQTL Summary\n")
  cat("---------------\n")
  cat("Tested probes    :",length(object),"\n")
  invisible(object)
} 
