## =============================================================================
## steady.2D -- solves the steady-state condition of
## ordinary differential equation systems
## 2-D reaction-transport models
## steady.band solves single-component 1-D reaction-transport models
## has similar calling sequence as integration routines from package deSolve
## =============================================================================

steady.2D    <- function (y,
                       time=0,
                       func,
                       parms=NULL,
                       nspec=NULL,
                       dimens,
                       ...)
{
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
    stop ("cannot run steady.2D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens)) 
    stop ("cannot run steady.2D: dimens should be specified")
  if (length(dimens)!=2) 
    stop ("cannot run steady.2D: dimens should contain 2 values")
  N     <- length(y)
  if (N%%prod(dimens) != 0) 
    stop("cannot run steady.2D: dimensions are not an integer fraction of number of state variables")
  if (is.null(nspec)) 
    nspec = N/prod(dimens)
  else if (nspec * prod(dimens) != N) 
    stop("cannot run steady.2D: dimens[1]*dimens[2]*nspec is not equal to number of state variables")
  

  out <- stodes(y=y,time=time,func=func,parms=parms,
                nnz=c(nspec,dimens),sparsetype="2D",...)                    
  return(out)
}

