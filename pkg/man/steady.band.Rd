\name{steady.band}
\alias{steady.band}
\title{Steady-state solver for ordinary differential equations; assumes a banded jacobian}
\description{Estimates the steady-state condition for a system of ordinary differential equations. 
  \cr Assumes a banded jacobian matrix, but does not rearrange the state variables. This is in contrast to \code{steady.1D}. 
Suitable for 1-D models that include transport only between adjacent layers and that model only one species}
\usage{steady.band(y, time=0, func, nspec=NULL, 
   bandup=nspec, banddown=nspec, ...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{time }{time for which steady-state is wanted; the default is \code{time}=0}
  \item{func }{either a user-supplied function that computes the values of the derivatives in the ode system (the model defininition) at time \code{time}, or a character string giving the name of a compiled function in a dynamically loaded shared library. 
     		See details}
  \item{nspec }{the number of *species* (components) in the model. }
  \item{bandup }{the number of nonzero bands above the jacobian diagonal}
  \item{banddown }{the number of nonzero bands below the jacobian diagonal}
  \item{... }{additional arguments passed to function \code{stode}}
}
\value{
  		A vector with the state variable values from the last iteration during estimation of steady-state condition of the system of equations, 
    plus the number of "global" values returned 
    
    If \code{y} has a names attribute, it will be used to label the output values. 
    The output will have the attribute \code{steady}, which returns TRUE, if steady-state has been reached and the attribute
    \code{precis} with the precision attained during each iteration.

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# 1000 simultaneous equations, solved 6 times for different 
# values of parameter "decay"
  model <- function (time,OC,parms,decay)
  {
   # model of particles (OC) that sink out of the water and decay
   # exponentially declining sinking rate, maximal 100 m/day
   sink <- 100*exp(-0.01*dist)
  
   # boundary flux; upper boundary=imposed concentration (100)
   Flux <- sink * c(100 ,OC)     
   
   # Rate of change= Flux gradient and first-order consumption
   dOC  <- -diff(Flux)/dx - decay*OC
   list(dOC,maxC=max(OC))
  }
  dx   <- 1                          # thickness of boxes
  dist <- seq(0,1000,by=dx)          # water depth at each modeled box interface

  ss   <- NULL 
  for (decay in seq(from=0.1,to=1.1,by=0.2))
    ss   <- cbind(ss,steady.band(runif(1000),func=model,
    parms=NULL,nspec=1,decay=decay)$y)  

  matplot(ss,1:1000,type="l",lwd=2,main="steady.band", ylim=c(1000,0),
          ylab="water depth, m", xlab="concentration of sinking particles")
  legend("bottomright",legend=seq(from=0.1,to=1.1,by=0.2),lty=1:10,
          title="decay rate",col=1:10,lwd=2)}
\details{This is the method of choice for single-species 1-D models. 

For multi-species 1-D models, this method can only be used if the state variables are arranged per box, per species (e.g. A[1],B[1],A[2],B[2],A[3],B[3],.... for species A, B). 

Usually a 1-D *model* function will have the species arranged as A[1],A[2],A[3],....B[1],B[2],B[3],.... 
in this case, use \code{steady.1D}}
\seealso{\code{\link{stode}} for the additional options 
   \item \code{\link{steady}}, for solving steady-state when the jacobian matrix is full 
  \item \code{\link{steady.1D}}, for solving steady-state, when the jacobian matrix is banded, and where the state variables need to be rearranged}
\keyword{math}

