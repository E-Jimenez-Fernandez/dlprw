\name{dlp_rw}
\alias{dlp_rw}
\title{Compute redundancy-weighted composite indicators}
\usage{
dlp_rw(x, polarity = NULL, p = 2, use_normat = TRUE,
       custom_min = NULL, custom_max = NULL,
       rank_tol = 1e-10, max_regressions = 1e6)
}
\arguments{
  \item{x}{Numeric matrix or data frame with observations in rows and nonconstant indicators in columns. Missing and nonfinite values are not accepted.}
  \item{polarity}{Names or indices of positive-polarity indicators; NULL means all positive. Other indicators are reversed.}
  \item{p}{One or several finite norm parameters greater than or equal to one.}
  \item{use_normat}{Logical. If FALSE, the input must already be oriented and scaled to [0,1], and polarity and custom bounds must be omitted.}
  \item{custom_min, custom_max}{Optional fixed normalization bounds, one per indicator.}
  \item{rank_tol}{Relative singular-value tolerance for numerical rank.}
  \item{max_regressions}{Upper bound on the number of exact nonempty subset calculations; Inf disables the safeguard.}
}
\description{Calculates normalized permutation-averaged redundancy weights and weighted lp composite scores. The weight calculation is independent of p.}
\details{The correction for an indicator is the average residual variance proportion over all subsets of preceding indicators, weighted by their permutation probabilities. The procedure measures linear redundancy rather than normative importance. See \code{\link{dlp_rw_weights}}.}
\value{A list with normalized data (\code{normat}), weights (\code{Weights}), correction factors (\code{q}), average linear redundancy (\code{redundancy}), weighted scores (\code{Comp.Indicator}), equally weighted scores (\code{Equal.Indicator}), descending ranks (\code{Ranking}, average ties), norm parameters (\code{p}), sample correlation matrix (\code{correlation}), normalization specification (\code{normalization}), number of subset calculations (\code{n_regressions}), rank tolerance (\code{rank_tol}) and method label (\code{method}).}
\examples{
set.seed(123)
x <- runif(50)
X <- data.frame(A = x, B = 0.7*x + 0.3*runif(50), C = runif(50))
fit <- dlp_rw(X, p = c(1, 2, 5))
fit$Weights
head(fit$Comp.Indicator)
}
\seealso{\code{\link{dlp_rw_weights}}, \code{\link{dlp_rw_scores}}, \code{\link{dlp_rw_normalize}}}
