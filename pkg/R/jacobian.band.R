
## =============================================================================
## jacobian.band  : multidiagonal (banded) jacobian matrix by differencing
## =============================================================================

jacobian.band<- function(y,         # (state) variables
                     func,          # function that calculates rate of change
                     bandup=1,      # number of nonzero bands above
                     banddown=1,    # number of nonzero bands below diagonal
                     dy=NULL,       # reference rate of change 
                     time=0,        # time passed to function 'func'
                     parms=NULL,    # parameter values passed to function 'func'
                     pert=1e-8,
                     ...)           # other arguments passed to function 'func'

{
# Reference value of state

  if (!is.numeric(y)) stop("y-values should be numeric")
  ny <- length(y)
  if(is.null(dy)) dy <-  unlist( func(time,y,parms,...))[1:ny]
  if (! is.numeric(dy)) stop("dy-values should either be NULL or numeric")
  if (length(dy) != ny) stop("function should return at least one value for each y")
   
  ynames  <- attr(y,"names")

  refy  <- y
  refdy <- dy

# Perturb the state variables, assume banded structure (only 3*nspec iterations)

  nband  <- bandup+banddown+1
  jacob  <- matrix(nrow=nband,ncol=ny,data=0)
  delt   <- perturb(y,pert)        # perturbation factors
  for ( j in 1:nband) {
    kpert    <- seq(j,ny,nband)                 # list of state var to perturb
    y[kpert] <- y[kpert] + delt[kpert]   # perturbed state var

    # new rate of change
    dy <-  unlist( func(time,y,parms,...))[1:ny]
    for (k in kpert) {
      iseq <- seq(max(k-bandup,1),min(k+banddown,ny))
      # impact of the selected states on the rate of change of all states
      jacob [iseq-k+bandup+1,k] <- (dy[iseq] -refdy[iseq]      )/delt[k]
    }
    y<-refy
  }
  colnames (jacob) <- ynames
  return(jacob )   # jacobian matrix, banded format

} ## END jacobian.band

