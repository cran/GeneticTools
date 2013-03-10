# Version: 30-11-2012, Daniel Fischer

summary.eqtl <- function(object, ...){
  
  cat("EQTL Summary - Testing, NO VALID OUTPUT!!!\n")
  cat("---------------\n")
  cat("Type of test     :",object$method,"\n")
  cat("Tested genes     :",length(object),"\n")
  invisible(object)
} 
