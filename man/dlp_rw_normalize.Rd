\name{dlp_rw_normalize}
\alias{dlp_rw_normalize}
\title{Orient and normalize indicators}
\usage{
dlp_rw_normalize(x, polarity = NULL, custom_min = NULL,
                 custom_max = NULL)
}
\arguments{
  \item{x}{Numeric matrix or data frame with finite observations.}
  \item{polarity}{Indices or names of positive-polarity indicators. NULL means all positive; integer(0) means all negative.}
  \item{custom_min, custom_max}{Optional numeric vectors of fixed bounds, aligned by names if supplied.}
}
\description{Applies min-max normalization, reversing negative-polarity attributes. Bounds must be strictly increasing and cover all observations; no clipping is performed.}
\value{A list containing normalized matrix Z, lower and upper bounds (min and max), and names of positive-polarity indicators (polarity).}
\examples{
dlp_rw_normalize(cbind(A = c(2, 4), B = c(10, 20)), polarity = "A")
}
