# DLp-RW: distance-based lp aggregation with permutation-averaged
# redundancy weights. Base R implementation of the manuscript's equations.

.dlp_matrix <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1))))
      stop("All indicators must be numeric.")
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L)
    stop("x must be a nonempty numeric matrix or data.frame.")
  if (any(!is.finite(x))) stop("Missing or nonfinite observations are not allowed.")
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  if (anyNA(colnames(x)) || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x)))
    stop("Indicator names must be nonempty and unique.")
  storage.mode(x) <- "double"
  x
}

.dlp_align <- function(x, nm, what) {
  if (!is.numeric(x) || length(x) != length(nm) || any(!is.finite(x)))
    stop(what, " must be a finite numeric vector of length m.")
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), nm))
      stop(what, " names must match the indicator names.")
    x <- x[nm]
  }
  setNames(as.numeric(x), nm)
}

