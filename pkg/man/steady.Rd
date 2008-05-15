\name{steady}
\alias{steady}
\title{General steady-state solver for a set of ordinary differential equations}
\description{Estimates the steady-state condition for a system of ordinary differential equations. 
	\cr This is a wrapper around steady-state solvers \code{stode} and \code{stodes}.}
\usage{steady(y, time=0, func, method="stode",...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{time }{time for which steady-state is wanted; the default is time=0}
  \item{func }{either a user-supplied function that computes the values of the derivatives in the ode system (the model defininition) at time \code{time}, 
    or a character string giving the name of a compiled function in a dynamically loaded shared library. 
     		See details}
  \item{method }{the solution method to use, one of \code{stode} or \code{stodes}}
  \item{... }{additional arguments passed to function \code{stode} or \code{stodes}}
}
\value{
  A vector with the state variable values from the last iteration during estimation of steady-state condition of the system of equations, plus the number of "global" values returned 
    
    If \code{y} has a names attribute, it will be used to label the output values. 
    The output will have the attribute \bold{steady}, which returns TRUE, if steady-state has been reached and the attribute
    \bold{precis} with the precision attained during each iteration.

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\seealso{\code{\link{stode}} and \code{\link{stodes}} for the additional options 
  \item \code{\link{steady.1D}}, \code{\link{steady.band}}, for solving steady-state when the jacobian matrix is banded}
\keyword{math}

