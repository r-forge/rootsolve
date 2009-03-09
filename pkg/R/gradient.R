
## =============================================================================
##                                                                            ##
##   routines that find the root of a nonlinear function                      ##
##                                                                            ##
## jacobian.full : generates a full jacobian matrix by numerical differencing ##
## jacobian.band : multidiagonal (banded) jacobian matrix by differencing     ##
## multiroot  : root of multiple simultaneous nonlinear equations             ##
##                                                                            ##
## internal functions - not to be included in package:                        ##
## -----------------------------                                              ##
## Perturb     : adds the numerical differencing value to a value             ##
##                                                                            ##
## =============================================================================

## =============================================================================
## gradient  : generates a full jacobian matrix by numerical differencing
## =============================================================================

gradient<- function(f,    # function returning a set of function values, as a vector
                    x,    # variables
                    centered = FALSE,
                    pert = 1e-8,
                    ...)  # additional arguments passed to function "f"...)              
{

## Reference value of variables and function value
  if (!is.numeric(x)) stop("x-values should be numeric")
       
  refx <- x
  reff <- f(x,...)

  Nx <- length(x)
  Nf <- length(reff)

## Perturb the state variables one by one
  delt   <- perturb(x,pert)
  jacob  <- matrix(nrow=Nf,ncol=Nx,data=0)

  for (j in 1:Nx) {
  # forward
    x[j] <- x[j]+delt[j]

     # recalculate model function value
    newf  <- f(x,...)
    del   <- (newf-reff)/delt[j]
           
    if (centered) {
    # backward formula
      x[j] <- refx[j]-delt[j]
      # recalculate model function value
      newf  <- f(x,...)
      del   <- (del-(newf-reff)/delt[j])/2
    }

    # impact of the current variable on function values
    jacob [,j] <- del

    x[j] <- refx[j]   # restore
  } # end for
  colnames(jacob) <- names(x)
  rownames(jacob) <- attr(del,"names")
  return(jacob)

} ## END gradient
