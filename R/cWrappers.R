# Version: 30-11-2012, Daniel Fischer

mdr.C <- function(X, fold, status,t,cv,cvp,top,na,fix){
  .Call( "mdr", X, fold, status, t, cv, cvp,top,na,fix)
} 
