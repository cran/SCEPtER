\name{stdGrid}
\docType{data}
\alias{stdGrid}
\title{Standard grid for mass and radius estimation}
\description{
The standard grid for mass and radius determination through a maximum likelihood
technique. The grid is computed assuming the solar scaled mixing-length
value.  
}
\usage{
stdGrid
}
\format{A data frames containing 9 variables:
  \describe{
    \item{Teff}{The effective temperature of the star (in K).}
    \item{logg}{The logarithm of the surface gravity of the star (in cm s^-2).}
    \item{FeH}{The metallicity [Fe/H] of the star.}
    \item{Dni}{The large frequency separation for the star (in unit of solar
      large frequency separation).}
    \item{nimax}{The frequency of maximum oscillation power for the star
      (in unit of solar frequency of maximum oscillation power).}
    \item{M}{The mass of the star (in unit of solar mass).}
    \item{R}{The radius of the star (in unit of solar radius)}
    \item{logAge}{The logarithm of the age of the star (in yr).}
    \item{pcAge}{The ratio between the age of the star and the time needed for
      central hydrogen depletion for that star. }
  }
}
\source{G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2013). \emph{Uncertainties in grid-based estimates of stellar mass and
  radius. SCEPtER: Stellar CharactEristics Pisa
  Estimation gRid}. In preparation.
}

\keyword{datasets}
