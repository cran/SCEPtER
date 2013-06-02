\name{estimatePert}
\alias{estimatePert}
\title{Mass and radius estimates for a single object}
\description{
The function computes the grid-based estimates of the mass and radius for
one object. A Monte Carlo simulation is performed to estimate the
statistical uncertainty on the derived quantities.
 }
\usage{
    estimatePert(data, STAR, sigma, thr, sel, Nrun=10000)
     }
\arguments{
\item{data}{the matrix estimation grid. }
\item{STAR}{a vector of 5 elements, containing the observational data.}
\item{sigma}{a vector of 5 elements, containing the error on observational
  data. The error on 
  seismic parameters are expressed as percentage.}
\item{thr}{the threshold for the selection of cases over which the
  computation of likelihood is performed.}
\item{sel}{a vector of 5 elements, containing 1 or 0 for inclusion or
  exclusion of a observational data from the likelihood computation.}
\item{Nrun}{the number of replication in the Monte Carlo simulation.}
}
\value{

  \code{estimatePert} returns a data frame obtained with the \code{Nrun}
  estimates
  of stellar parameters. For each simulation, the object to be
  reconstructed is obtained sampling from a
  multivariate normal distribution with mean given by \code{STAR} and
  covariance matrix \code{diag(sigma^2)}.

  Mass and radius estimates are obtained by a maximum likelihood
  technique widely adopted in literature. Details on the technique can
  be found in the references reported below.
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
# Solar seismic data from  Thiery et al. (2000) 
Dnisun <- 134.8
nimaxsun <- 3034

# Observational data for alpha Cent B from:
# Kjeldsen et al. (2005), Porto de Mello et al. (2008)
# Observed mass and radius: 0.935 +- 0.006, 0.863 +- 0.005
# from Miglio & Montalban (2005) 
acenB <- c(5316, 0, 0.25, 161.5/Dnisun, 4100/nimaxsun)
sigmaCen <- c(30, 0, 0.04, 0.01, 0.05)

# Load the standard estimation grid
data(stdGrid)

# For accurate estimates adopt at least Nrun=10000
res <- estimatePert(stdGrid, acenB, sigmaCen, 3, c(1,0,1,1,1), Nrun=200)

quantile(res$M, c(0.16, 0.5, 0.84))
quantile(res$R, c(0.16, 0.5, 0.84))

}
\keyword{manip}
