\name{errorObs}
\alias{errorObs}
\title{Perturbe a set of observational data}
\description{
The function returns a perturbed object, starting from observational data
and a vector of uncertainties. 
 }
\usage{
    errorObs(sigma, STAR)
     }
\arguments{
\item{sigma}{a vector of 5 elements containing the uncertainties on
  observational data. The error on 
  seismic parameters must be expressed as percentage.} 
\item{STAR}{a vector of 5 elements containing  the observational data.}
}
\value{
  \code{errorObs} returns an object obtained sampling from a
  multivariate normal distribution with vector of mean \code{star} and
  covariance matrix \code{diag(sigma^2)}. 
}
\details{
  The vector of observations \code{STAR} must contain, in the order:
  \enumerate{
  \item{The effective temperature (in K) of the star.}
  \item{The logarithm of the surface gravity of the star (in cm s^-2).}
  \item{The metallicity [Fe/H] of the star.}
  \item{The value of large frequency separation of the star, divided by
    the solar value of this quantity. }
  \item{The value of the frequency of maximum oscillation power of the
  star, divided by 
    the solar value of this quantity. }
}

  The vector of the uncertainties on the observation \code{sigma} must
  contain the uncertainty on the above quantity. The
  uncertainty on the last 
  two values must be expressed as percentage.
  }
\keyword{manip}
