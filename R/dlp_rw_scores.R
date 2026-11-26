#' Weighted lp scores for fixed weights (one or several values of p)
#' @param Z Normalized data, with columns matching the weights.
#' @param weights Strictly positive weights summing to one.
#' @param p Finite norm parameter(s), each >= 1.
#' @return Matrix with one column per p, in input order.
dlp_rw_scores <- function(Z, weights, p = 2) {
  Z <- .dlp_matrix(Z)
  if (any(Z < 0 | Z > 1)) stop("Z must belong to [0,1].")
  w <- .dlp_align(weights, colnames(Z), "weights")
  if (any(w <= 0) || abs(sum(w) - 1) > 1e-10)
    stop("Weights must be strictly positive and sum to one.")
  w <- w / sum(w) # correct floating-point summation only
  if (!is.numeric(p) || !length(p) || any(!is.finite(p) | p < 1))
    stop("p must contain finite numbers >= 1.")
  # Row scaling avoids underflow for large p and small nonzero observations.
  a <- apply(Z, 1, max)
  B <- sweep(Z, 1, ifelse(a > 0, a, 1), "/")
  out <- vapply(p, function(power) a * drop((B^power) %*% w)^(1/power),
                numeric(nrow(Z)))
  out <- matrix(out, nrow = nrow(Z), ncol = length(p),
                dimnames = list(rownames(Z), paste0("p_", format(p, trim = TRUE))))
  out
}

