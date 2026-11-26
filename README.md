# dlprw

Permutation-averaged linear redundancy weighting for norm-based composite indicators (DLp-RW).

The method averages the Pena-Trapero correction factors over all orders of the indicators, using equivalent subset enumeration. It assigns positive normalized weights independently of the aggregation parameter `p`. The method corrects linear statistical redundancy, not normative importance; it does not guarantee invariance under arbitrary duplications or detect every nonlinear dependence.

## Installation

```r
install.packages("remotes")
remotes::install_github("E-Jimenez-Fernandez/dlprw")
library(dlprw)
```

This is a development package (version 0.1.0). The repository must be published before the GitHub command will work.

## Example

```r
set.seed(123)
n <- 100
income <- runif(n)
X <- data.frame(
  Income = income,
  Education = 0.7 * income + 0.3 * runif(n),
  Health = runif(n),
  Employment = runif(n)
)
fit <- dlp_rw(X, polarity = 1:4, p = c(1, 2, 3, 5, 10))
fit$Weights
fit$q
head(fit$Comp.Indicator)
head(fit$Ranking)
cor(fit$Comp.Indicator, method = "kendall")
```

All indicators must be numeric, nonconstant and finite. `polarity` identifies positive indicators; others are reversed. If omitted, all are positive. For already normalized and oriented data, set `use_normat = FALSE`. Exact subset enumeration scales exponentially in the number of attributes; `max_regressions` limits its cost. SVD uses a relative numerical rank tolerance (`rank_tol`).

## Functions

- `dlp_rw()` — complete procedure.
- `dlp_rw_normalize()` — min-max normalization and polarity.
- `dlp_rw_weights()` — exact subset-averaged redundancy weights.
- `dlp_rw_scores()` — scores for fixed weights and one or more p values.

Use `?dlp_rw` and the help pages for details. No third-party R packages are required for computation.

## Reproduction

The `reproduction/` directory contains the R and Python Monte Carlo scripts used in the manuscript development. Install this R package before running the R script. See the directory README for execution notes. Python and R random draws are not identical.

## Development

```r
install.packages("devtools")
devtools::check()
```

The package tests cover positivity, normalization, equivariance, duplication cases, and agreement with explicit permutation enumeration.

## License

MIT (see LICENSE).
