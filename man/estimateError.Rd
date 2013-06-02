\name{estimateError}
\alias{estimateError}
\title{Monte Carlo simulation for a single object}
\description{
The function computes the grid-based estimates of the mass and radius for
several objects perturbing their observational data with a random
uncertainty.
It can be used to assess the performance of the recovery procedure
sampling from different grids.
 }
\usage{
    estimateError(data, STAR, sigma, thr, sel)
     }
\arguments{
  \item{data}{the matrix estimation grid. }
\item{STAR}{a matrix of 5 columns, containing by row the observational
  data of the objects.}
\item{sigma}{a vector of 5 elements, containing the error on observational
  data. The error on 
  seismic parameters are expressed as percentage.}
\item{thr}{the threshold for the selection of cases over which the
  computation of likelihood is performed}
\item{sel}{a vector of 5 elements, containing 1 or 0 for inclusion or
  exclusion of a observational quantity from the likelihood computation.}
}
\value{

  \code{estimateError} returns a data frame with the
  grid-based estimates
  of stellar parameters for each row in the matrix \code{STAR}.
  In the columns of the output object there are:
  \enumerate{
    \item{The grid-based estimate of the mass of the star.}
    \item{The grid-based estimate of the radius of the star.}
    \item{The number of points in the sampling grid which concurred in
      the estimates of mass and radius. For details on
      the technique, see references below.}
    \item{ The true mass of the star.}
    \item{ The true radius of the star.}
    \item{The logarithm of the age of the star (in yr).}
    \item{The ratio between the age of the star and the time needed for
      central hydrogen depletion for that star.}
    \item{The metallicity [Fe/H] of the star.}
    \item{The relative error on the grid-based estimate of the mass.}
    \item{The relative error on the grid-based estimate of the radius.}
    }
  
  For each row \eqn{i} in \code{STAR}, the observational data are
  perturbed adding a random noise, sampled from a 
  multivariate normal distribution with mean given by the \eqn{i}-th
  row of \code{STAR} and 
  covariance matrix \code{diag(sigma^2)}. 
}
\references{
  
G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2013). \emph{Uncertainties in grid-based estimates of stellar mass and
  radius. SCEPtER: Stellar CharactEristics Pisa
  Estimation gRid}. In preparation.

S. Basu, G.A. Verner, W.J. Chaplin, and Y. Elsworth (2012). \emph{Effect
  of Uncertainties in Stellar Model Parameters on Estimated Masses and
  Radii of Single Stars}.  The
Astrophysical Journal, 746, 76.

N. Gai, S. Basu, W.J. Chaplin, and Y. Elsworth (2011). \emph{An In-depth
  Study of Grid-based Asteroseismic Analysis}.  The
Astrophysical Journal, 730, 63.

}

\examples{
data(stdGrid)

sigma <- c(100, 0.25, 0.1, 0.025, 0.05)
star <- sampleStar(200, stdGrid)

res <- estimateError(stdGrid, star, sigma, 3, c(1, 0, 1, 1, 1))

plot(density(res$errorM))

plot(R ~ trueR, data=res)
}
\keyword{manip}
