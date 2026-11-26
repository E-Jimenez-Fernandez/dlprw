#' Exact permutation-averaged redundancy weights
#' @param Z Normalized numeric matrix with nonconstant columns.
#' @param rank_tol Relative singular-value tolerance for numerical rank.
#' @param max_regressions Computational safeguard (Inf disables it).
#' @return Named q factors and weights, computation diagnostics.
dlp_rw_weights <- function(Z, rank_tol = 1e-10, max_regressions = 1e6) {
  Z <- .dlp_matrix(Z)
  n <- nrow(Z); m <- ncol(Z); nm <- colnames(Z)
  if (n < 2L) stop("At least two observations are required.")
  if (any(Z < 0 | Z > 1)) stop("Z must belong to [0,1].")
  if (!is.numeric(rank_tol) || length(rank_tol) != 1L ||
      !is.finite(rank_tol) || rank_tol <= 0 || rank_tol >= 1)
    stop("rank_tol must belong to (0,1).")
  if (!is.numeric(max_regressions) || length(max_regressions) != 1L ||
      is.na(max_regressions) || max_regressions < 0)
    stop("max_regressions must be nonnegative.")
  count <- m * (2^(m - 1) - 1)
  if (!is.finite(count) || count > max_regressions)
    stop("Exact subset calculation exceeds max_regressions: ", count,
         ". Reduce m or explicitly increase the computational limit.")

  # Center each indicator and normalize its Euclidean length.
  # Column scaling does not change the theoretical predictor subspaces.
  V <- sweep(Z, 2, colMeans(Z), "-")
  # Pre-scale by maximum absolute deviation to avoid underflow.
  scale <- apply(abs(V), 2, max)
  # Test exact nonconstancy on original columns as well.
  constant <- apply(Z, 2, function(z) all(z == z[1L])) | scale == 0
  if (any(constant)) stop("Constant columns: ", paste(nm[constant], collapse = ", "))
  V <- sweep(V, 2, scale, "/")
  U <- sweep(V, 2, sqrt(colSums(V^2)), "/")

  # For each j, the empty subset contributes exactly 1/m.
  q <- setNames(rep(1/m, m), nm)
  if (m > 1L) {
    for (j in seq_len(m)) {
      others <- setdiff(seq_len(m), j)
      y <- U[, j]
      for (k in seq_len(m - 1L)) {
        # Each subset of size k receives probability 1/(m*choose(m-1,k)).
        alpha <- 1 / (m * choose(m - 1L, k))
        # Compute residual norm using the SVD projection, not matrix inversion.
        contributions <- utils::combn(seq_along(others), k, FUN = function(idx) {
          S <- others[idx]
          sv <- base::svd(U[, S, drop = FALSE], nv = 0)
          keep <- which(sv$d > rank_tol * sv$d[1L])
          if (length(keep)) {
            basis <- sv$u[, keep, drop = FALSE]
            residual <- y - drop(basis %*% crossprod(basis, y))
          } else residual <- y
          # 1-R^2 = ||(I-P_S)y||^2 because ||y|| = 1.
          min(1, max(0, sum(residual^2)))
        })
        q[j] <- q[j] + alpha * sum(contributions)
      }
    }
  }
  w <- q / sum(q)
  list(q = q, weights = w, n_regressions = count,
       rank_tol = rank_tol, method = "exact subset averaging")
}

