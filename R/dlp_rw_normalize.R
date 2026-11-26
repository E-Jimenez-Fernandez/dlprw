#' Polarity-based min-max normalization
#' @param x Numeric matrix or data frame.
#' @param polarity Indices or names of positive-polarity indicators.
#' NULL means all positive; integer(0) means all negative.
#' @param custom_min,custom_max Optional fixed normalization bounds.
#' @return A normalized matrix and the normalization specification.
dlp_rw_normalize <- function(x, polarity = NULL, custom_min = NULL,
                             custom_max = NULL) {
  X <- .dlp_matrix(x)
  nm <- colnames(X); m <- ncol(X)
  if (is.null(polarity)) polarity <- seq_len(m)
  if (is.character(polarity)) {
    if (anyNA(polarity) || !all(polarity %in% nm)) stop("Unknown polarity names.")
    polarity <- match(polarity, nm)
  }
  if (!is.numeric(polarity) || anyNA(polarity) ||
      any(!is.finite(polarity)) || any(polarity != floor(polarity)) ||
      any(polarity < 1 | polarity > m)) stop("Invalid positive-polarity indices.")
  lo <- if (is.null(custom_min)) apply(X, 2, min) else .dlp_align(custom_min, nm, "custom_min")
  hi <- if (is.null(custom_max)) apply(X, 2, max) else .dlp_align(custom_max, nm, "custom_max")
  span <- hi - lo
  if (any(!is.finite(span) | span <= 0))
    stop("Normalization requires strictly increasing finite bounds; check constant columns.")
  if (any(sweep(X, 2, lo, "-") < 0) || any(sweep(X, 2, hi, "-") > 0))
    stop("Observations outside normalization bounds; no clipping is applied.")
  Z <- sweep(sweep(X, 2, lo, "-"), 2, span, "/")
  neg <- setdiff(seq_len(m), polarity)
  if (length(neg)) Z[, neg] <- 1 - Z[, neg, drop = FALSE]
  list(Z = Z, min = lo, max = hi, polarity = nm[sort(unique(polarity))])
}

