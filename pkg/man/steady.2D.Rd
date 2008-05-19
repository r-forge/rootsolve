\name{steady.2D}
\alias{steady.2D}
\title{Steady-state solver for multicomponent 1-D ordinary differential equations}
\description{Estimates the steady-state condition for a system of ordinary differential equations. 
  \cr Assumes a sparse jacobian matrix. 
  \cr Suitable for 2-D models that include transport only between adjacent layers }
\usage{steady.2D(y, time=0, func, dimens=NULL,...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{time }{time for which steady-state is wanted; the default is time=0}
  \item{func }{either a user-supplied function that computes the values of the derivatives in the ode system (the model defininition) at time \code{time}, 
    or a character string giving the name of a compiled function in a dynamically loaded shared library. 
     		See details}
  \item{dimens}{a 2-valued vector with the dimensionality of the model, i.e. the number of *boxes* in x- and y-direction}
  \item{... }{additional arguments passed to function \code{stodes}}
}
\value{
  A vector with the state variable values from the last iteration during estimation of steady-state condition of the system of equations, plus the number of "global" values returned 
    
    If \code{y} has a names attribute, it will be used to label the output values. 
    The output will have the attribute \bold{steady}, which returns TRUE, if steady-state has been reached and the attribute
    \bold{precis} with the precision attained during each iteration.

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
###########################################
# Diffusion in 2-D; imposed boundary conditions
diffusion2D <- function(t,Y,par)
  {
   y    <- matrix(nr=n,nc=n,data=Y)
   dY   <- -r*y        # consumption
   BND   <- rep(1,n)   # boundary concentration 

   #diffusion in X-direction; boundaries=imposed concentration
   Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
   dY   <- dY - (Flux[2:(n+1),]-Flux[1:n,])/dx

   #diffusion in Y-direction
   Flux <- -Dy * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dy
   dY    <- dY - (Flux[,2:(n+1)]-Flux[,1:n])/dy

   return(list(as.vector(dY)))
  }

  # parameters
  dy    <- dx <- 1   # grid size
  Dy    <- Dx <- 1   # diffusion coeff, X- and Y-direction
  r     <- 0.025     # consumption rate

  n  <- 100
  y  <- matrix(nr=n,nc=n,10.)

  ST3 <- steady.2D(y,func=diffusion2D,parms=NULL,pos=TRUE,dimens=c(n,n),
                 lrw=1000000,atol=1e-10,rtol=1e-10,ctol=1e-10)
  y <- matrix(nr=n,nc=n,data=ST3$y)
    filled.contour(y,color.palette=terrain.colors)
    
}
\details{}
\seealso{\code{\link{stodes}} for the additional options 
  \item \code{\link{steady}}, for solving steady-state when the jacobian matrix is full 
  \item \code{\link{steady.1D}}, for solving steady-state for 1-D models }
\keyword{math}

