% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirichlet.R
\name{rdirichlet}
\alias{rdirichlet}
\alias{rdirichlet1}
\alias{ddirichlet}
\title{Dirichlet distribution}
\usage{
rdirichlet(n, alpha)

rdirichlet1(alpha)

ddirichlet(x, alpha)
}
\arguments{
\item{n}{Number of random vectors to generate.}

\item{alpha}{Vector or (for ddirichlet) matrix containing shape parameters.}

\item{x}{A vector containing a single random deviate or matrix containg one
random deviate per row.}
}
\description{
Density and random generation for the Dirichlet distribution with mean equal
to mean and standard deviation equal to sd.
}
\details{
The Dirichlet is the multidimensional of the bet distribution and the
canonical Bayesian Distribution for the parameters of a multinomial
distribution. \code{rdirichlet1(alpha)} is a shortcut for
\code{rdirichlet(1, alpha)}
}
\examples{

set.seed(42)
rdirichlet(5, c(1,1,1) )
set.seed(42)
rdirichlet1(c(1,1,1))

}
\seealso{
\code{\link{dbeta}}, \code{\link{rbeta}}
}
\author{
Code original posted by Ben Bolker to R-News on Fri Dec 15 2000. See
https://stat.ethz.ch/pipermail/r-help/2000-December/009561.html. Ben
attributed the code to Ian Wilson i.wilson@maths.abdn.ac.uk. Subsequent
modifications by Gregory R. Warnes greg@warnes.net.
}
\keyword{cluster}
\keyword{~kwd2}
