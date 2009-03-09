
## =============================================================================
## uniroot.all: multiple roots of one nonlinear equation
## =============================================================================

uniroot.all <- function(
        f,                    # the function for which the root is sought.
        interval,             # a vector containing the end-points of the interval to be searched for the root.
        lower= min(interval), # the lower end point of the interval to be searched.
        upper= max(interval), # the upper end point of the interval to be searched.
        tol= .Machine$double.eps^0.2,      # the desired accuracy (convergence tolerance).
        maxiter= 1000,       # the maximum number of iterations.
        n = 100,             #number of subintervals in which root is sought
        ... )               #additional named or unnamed arguments to be passed to f (but beware of partial matching to other arguments).
{
## error checking as in uniroot...
  if (!missing(interval) && length(interval) != 2)
     stop("'interval' must be a vector of length 2")
  if (!is.numeric(lower) || !is.numeric(upper) || lower >=
     upper)
    stop("lower < upper  is not fulfilled")
 
## subdivide interval in n subintervals and estimate the function values
  xseq <- seq(lower,upper,len=n+1)
  mod  <- f(xseq,...)

## some function values may already be 0
  Equi <- xseq[which(mod==0)]

  ss   <- mod[1:n]*mod[2:(n+1)]  # interval where functionvalues change sign
  ii   <- which(ss<0)

  for (i in ii)
    Equi <- c(Equi,uniroot(f,lower=xseq[i],upper=xseq[i+1],...)$root)

  return(Equi)
}

