\name{estimate}
\alias{estimate}
\title{Mass, radius, and age estimates for the supplied objects}
\description{
The function computes the grid-based estimates of the mass, radius, and age for
the supplied objects. 
 }
\usage{
    estimate(data, STAR, sigma, thr, sel, parallel=FALSE)
     }
\arguments{
\item{data}{the matrix estimation grid. The matrix should be sorted according to
  the values in the first column (the stellar effective
  temperature). }
\item{STAR}{a matrix of 7 columns, containing in each row the
  observational data of the stellar objects.} 
\item{sigma}{a vector of 7 elements, containing the error on observational
  data. The error on 
  seismic parameters are expressed as relative error.}
\item{thr}{the threshold for the selection of cases over which the
  computation of likelihood is performed.}
\item{sel}{a vector of 7 elements, containing 1 or 0 for inclusion or
  exclusion of a observational data from the likelihood computation.}
\item{parallel}{logical flag indicating if the estimations should be
  performed in parallel on a multi-core architecture.}
}
\value{
  \code{STAR} should contains in order: the stellar effective temperature,
  the logarithm of the surface gravity, the metallicity [Fe/H], the
  average large frequency separation (in micro Hz), the maximum
  oscillation power (in micro Hz), the mass (in solar units), the radius
  (in solar units). \code{sigma} should contain, in the same order, the
  uncertainties on the observational quantities.  
  
  \code{estimate} returns a data frame obtained with the
  estimates
  of stellar parameters for each row of the input matrix \code{STAR}.
  In the columns of the data frame there are: the identifier \code{id}
  of the row of the input matrix adopted for the estimate; the mass
  estimate \code{M} (in solar units); the radius estimate \code{R} (in
  solar units); the age estimate \code{age} (in Gyr).

  Mass, radius, and age estimates are obtained by a maximum likelihood
  technique widely adopted in literature. Details on the technique can
  be found in the references reported below.
}
\references{
  
G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2014). \emph{Uncertainties in grid-based estimates of stellar mass and
  radius. SCEPtER: Stellar CharactEristics Pisa
  Estimation gRid}. A&A, 561, A125.

G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2014). \emph{Uncertainties in asteroseismic grid-based estimates of stellar ages. SCEPtER: Stellar CharactEristics Pisa
  Estimation gRid}. A&A, accepted.

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
acenB <- c(5316, 0, 0.25, 161.5/Dnisun, 4100/nimaxsun, 0.935, 0.863)
sigmaCen <- c(30, 0, 0.04, 0.01, 0.05, 0.006/0.935, 0.005/0.863)

# Load the standard estimation grid
data(stdGrid)

# order standard grid by effective temperature (required by estimate)
ii <- do.call(order, as.data.frame(stdGrid[,1]))
oGrid <- stdGrid[ii, ]

# For accurate estimates adopt at least NUM=10000
NUM <- 200
# obtain NUM copies of the object and add Gaussian perturbations
star <- matrix(acenB, byrow=TRUE, nrow=NUM, ncol=7)
starp <- errorObs(sigmaCen, star)

# select the observational constraints
# using Teff, [Fe/H], Delta nu, nu_max 
vsel <- c(1,0,1,1,1,0,0)

# call the estimation routines
res <- estimate(oGrid, starp, sigmaCen, 3, vsel)

quantile(res$M, c(0.16, 0.5, 0.84))
quantile(res$R, c(0.16, 0.5, 0.84))

}
\keyword{manip}
