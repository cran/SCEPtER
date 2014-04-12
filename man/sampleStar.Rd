\name{sampleStar}
\alias{sampleStar}
\title{Sample a set of objects out of a grid}
\description{
The function extracts, without replacement, a sample of objects out of a grid.  
 }
\usage{
    sampleStar(n, grid, restrict=TRUE)
     }
\arguments{
\item{n}{the number of objects to be sampled.}
\item{grid}{the grid from which the objects are sampled.}
\item{restrict}{logical flag indicating if only models younger than 14 Gyr
  should be considered.}
}
\value{

  \code{sampleStar} returns a matrix of \code{n} rows, sampled without
  replacement from \code{grid}.

}

\examples{
data(stdGrid)
STARS <- sampleStar(10, stdGrid)
}
\keyword{manip}
