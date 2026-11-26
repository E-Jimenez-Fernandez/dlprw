\name{dlp_rw_scores}
\alias{dlp_rw_scores}
\title{Compute weighted lp composite scores}
\usage{
dlp_rw_scores(Z, weights, p = 2)
}
\arguments{
  \item{Z}{Numeric matrix or data frame with finite entries in [0,1].}
  \item{weights}{Strictly positive numeric vector summing to one, aligned with the columns of Z; named vectors are matched by name.}
  \item{p}{One or several finite norm parameters greater than or equal to one.}
}
\description{Computes weighted lp aggregation at fixed weights, allowing comparison of different aggregation geometries without re-estimating weights.}
\value{A numeric matrix of composite scores, one column per p.}
\examples{
Z <- rbind(c(0.5, 0.5), c(0.8, 0.2))
dlp_rw_scores(Z, c(0.5, 0.5), p = c(1, 2, 5))
}
