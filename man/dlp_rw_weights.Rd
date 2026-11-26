\name{dlp_rw_weights}
\alias{dlp_rw_weights}
\title{Compute exact permutation-averaged redundancy weights}
\usage{
dlp_rw_weights(Z, rank_tol = 1e-10, max_regressions = 1e6)
}
\arguments{
  \item{Z}{Numeric matrix or data frame in [0,1] with at least two rows and nonconstant columns.}
  \item{rank_tol}{Relative singular-value tolerance in (0,1).}
  \item{max_regressions}{Maximum number of nonempty subset calculations; Inf disables the safeguard.}
}
\description{Computes the average proportion of unexplained variance for each indicator across all possible orders, using centered orthogonal projections and singular value decomposition.}
\details{The empty subset contributes 1/m. Each subset of size k excluding the target indicator has probability 1/(m * choose(m-1,k)). The correction factors are normalized to sum to one. Exact enumeration is exponential in the number of indicators. Near rank deficiency the result may depend on the numerical rank tolerance. The method does not detect all nonlinear dependence.}
\value{A list containing named correction factors (q), named normalized weights (weights), number of subset calculations (n_regressions), rank_tol and method.}
\examples{
Z <- cbind(A = c(0,0,1,1), B = c(0,1,0,1))
dlp_rw_weights(Z)
dlp_rw_weights(cbind(Z, A_copy = Z[,1]))
}
