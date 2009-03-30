
## =============================================================================
## multiroot: root of multiple simultaneous nonlinear equations
## =============================================================================

multiroot <- function(f, start, maxiter=100, rtol=1e-6, atol=1e-8, ctol=1e-8,        # minimal change in dy
                      useFortran=TRUE, positive=FALSE, ...)  {

  N        <- length(start)
  if (!is.numeric(start))
    stop("start conditions should be numeric")
  if (!is.numeric(maxiter))
    stop("`maxiter' must be numeric")
  if (as.integer(maxiter) < 1)
    stop ("maxiter must be >=1")
  if (!is.numeric(rtol))
    stop("`rtol' must be numeric")
  if (!is.numeric(atol))
    stop("`atol' must be numeric")
  if (!is.numeric(ctol))
    stop("`ctol' must be numeric")
  if (length(atol) > 1 && length(atol) != N)
    stop("`atol' must either be a scalar, or as long as `start'")
  if (length(rtol) > 1 && length(rtol) != N)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (length(ctol) > 1)
    stop("`ctol' must be a scalar")

  if (useFortran) {
    Fun <- function (time=0,x,parms=NULL)
      list(f(x,...))
 
    x <- steady(y=start,time=0,func=Fun,parms=NULL,atol=atol,positive=positive,
              rtol=rtol,ctol=ctol,jacfunc="fullint",maxiter=maxiter)
    precis <- attr(x,"precis")
    attributes(x)<-NULL

    x        <- unlist(x)
    reffx    <-f(x,...)
    i        <- length(precis)
    
  } else {  # simple R-implementation

    precis   <- NULL

    x        <- start
    jacob    <- matrix(nrow=N,ncol=N,data=0)
    reffx    <- f(x,...)     # function value,
  
    if (length (reffx) != N)
      stop("'f', function must return as many function values as elements in start")


    for (i in 1:maxiter) {
      refx   <- x
      pp     <- mean(abs(reffx)) # check convergence...
      precis <- c(precis,pp)
      ewt    <- rtol*abs(x)+atol
      if (max(abs(reffx/ewt))<1) break

      # estimate jacobian
      delt   <- perturb(x)

      for (j in 1:N) {
        x[j] <- x[j]+delt[j]   # Perturb the state variables one by one
        fx   <-  f(x,...)      # recalculate function value

        # impact of the current state var on rate of change of all state vars
        jacob [,j] <- (fx-reffx)/delt[j]

        x[j] <- refx[j]   # restore
      }

   # new estimate 
      relchange <- as.numeric(solve(jacob,-1*reffx))
      if (max(abs(relchange)) < ctol) break
      x  <- x + relchange
      reffx  <- f(x,...)     # function value,
    
    } # end for
  } # end fortran/R
  names(x) <- names(start)

  return(list(root=x,f.root=reffx,iter=i,estim.precis=precis[length(precis)]))
}  # end multiroot

