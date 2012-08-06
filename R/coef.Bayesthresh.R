# S3 method for coeficients

coef.Bayesthresh <- function(object, ...)
  {
     if(!inherits(object, "Bayesthresh"))
       stop("Use an object of class Bayesthresh")
     fixed <- object$EfFixef
     return(fixed)
   }


