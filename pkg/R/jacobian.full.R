
## =============================================================================
## jacobian.full  : generates a full jacobian matrix by numerical differencing
## =============================================================================

jacobian.full<- function(y,             # (state) variables
                         func,          # function that calculates rate of change
                         dy=NULL,       # reference rate of change
                         time=0,        # time passed to function 'func'
                         parms=NULL,    # parameter values passed to function 'func'
                         pert=1e-8,
                         ...)           # other arguments passed to function 'func'
{
# Reference value of state variable and of rate of change (dC)
  if (!is.numeric(y)) stop("y-values should be numeric")
  N    <- length(y)
  refy <- y
  ifelse (is.null(dy), refdy <- unlist( func(time,y,parms,...))[1:N], refdy <- dy)
  if (! is.numeric(refdy)) stop("dy-values should either be NULL or numeric")
  if (length(refdy) != N)  stop("function should return at least one value for each y")
  ynames  <- attr(y,"names")
 
# Perturb the state variables one by one
  delt   <- perturb(y,pert)

  ny <-length(y)
  jacob  <- matrix(nrow=ny,ncol=ny,data=0)
  for (j in 1:ny) {
    y[j] <- y[j]+delt[j]

    # recalculate model rate of change
    dy  <-  unlist( func(time,y,parms,...))[1:N]

    # impact of the current state var on rate of change of all state vars
    jacob [,j] <- (dy-refdy)/delt[j]

    y[j] <- refy[j]   # restore
  }
  colnames (jacob) <- ynames
           
  return(jacob)

} ## END jacobian.full
