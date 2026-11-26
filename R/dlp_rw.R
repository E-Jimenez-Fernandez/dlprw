#' Main entry point for the new methodology
#' @param x Raw numeric indicators (or normalized data if use_normat=FALSE).
#' @param polarity Positive-polarity names/indices; NULL means all positive.
#' @param p One or several finite norm parameters >=1.
#' @param use_normat Apply polarity-based min-max normalization.
#' @param custom_min,custom_max Fixed normalization bounds, if desired.
#' @param rank_tol Relative numerical-rank tolerance for SVD.
#' @param max_regressions Maximum exact subset regressions.
#' @return A list with weights, redundancy factors, scores and rankings.
#' @export
dlp_rw <- function(x, polarity = NULL, p = 2, use_normat = TRUE,
                   custom_min = NULL, custom_max = NULL,
                   rank_tol = 1e-10, max_regressions = 1e6) {
  if (!is.logical(use_normat) || length(use_normat) != 1L || is.na(use_normat))
    stop("use_normat must be TRUE or FALSE.")
  if (use_normat) {
    norm <- dlp_rw_normalize(x, polarity, custom_min, custom_max)
    Z <- norm$Z
    spec <- norm[c("min", "max", "polarity")]
  } else {
    if (!is.null(polarity) || !is.null(custom_min) || !is.null(custom_max))
      stop("With use_normat=FALSE, input must already be oriented and normalized.")
    Z <- .dlp_matrix(x); spec <- NULL
  }
  fit <- dlp_rw_weights(Z, rank_tol, max_regressions)
  scores <- dlp_rw_scores(Z, fit$weights, p)
  equal <- dlp_rw_scores(Z, rep(1/ncol(Z), ncol(Z)), p)
  ranks <- vapply(seq_len(ncol(scores)), function(j) rank(-scores[, j], ties.method = "average"),
                  numeric(nrow(Z)))
  ranks <- matrix(ranks, nrow = nrow(Z), dimnames = dimnames(scores))
  list(normat = Z, Weights = fit$weights, q = fit$q,
       redundancy = 1 - fit$q, Comp.Indicator = scores,
       Equal.Indicator = equal, Ranking = ranks, p = p,
       correlation = stats::cor(Z), normalization = spec,
       n_regressions = fit$n_regressions, rank_tol = fit$rank_tol,
       method = fit$method)
}
