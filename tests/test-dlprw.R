library(dlprw)
near <- function(a, b, tol = 1e-8) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
# Balanced, orthogonal columns with positive variance.
Z <- cbind(A = c(0,0,1,1), B = c(0,1,0,1))
fit <- dlp_rw(Z, p = c(1,2,5), use_normat = FALSE)
stopifnot(near(fit$Weights, c(.5,.5)), near(fit$q, c(1,1)))
# Exact duplicate and unrelated indicator: score invariance.
D <- cbind(Z, A_copy = Z[, "A"])
dup <- dlp_rw(D, p = c(1,2,5), use_normat = FALSE)
stopifnot(near(dup$Weights, c(.25,.5,.25)), near(dup$Comp.Indicator, fit$Comp.Indicator))
# Permutation equivariance.
revfit <- dlp_rw(D[, c(3,1,2)], p = c(1,2,5), use_normat = FALSE)
stopifnot(near(revfit$Weights[names(dup$Weights)], dup$Weights),
          near(revfit$Comp.Indicator, dup$Comp.Indicator))
# All identical columns: equal weights, positive despite singularity.
ident <- dlp_rw(cbind(A=Z[,1], B=Z[,1], C=Z[,1]), use_normat=FALSE)
stopifnot(near(ident$q, rep(1/3,3)), near(ident$Weights, rep(1/3,3)))
# Correlated duplication: expected lack of exact invariance.
X <- cbind(A=c(0,0,1,1), B=c(0,0.3,0.5,1))
r2 <- cor(X)[1,2]^2
C <- dlp_rw(cbind(X, A_copy=X[,1]), use_normat=FALSE)
stopifnot(near(C$q, c((1-r2/3)/2,1-2*r2/3,(1-r2/3)/2)),
          sum(C$Weights[c(1,3)]) > .5)
# Positivity, normalization, number of regressions, p independence.
set.seed(42)
X <- matrix(runif(400), ncol=4, dimnames=list(NULL,LETTERS[1:4]))
f <- dlp_rw(X, p=c(1,2,5), use_normat=FALSE)
stopifnot(f$n_regressions==28, all(f$Weights >= 1/16-1e-10),
          near(sum(f$Weights),1),
          near(f$Weights, dlp_rw(X,p=3,use_normat=FALSE)$Weights),
          all(f$Comp.Indicator[,1] <= f$Comp.Indicator[,2]+1e-10),
          all(f$Comp.Indicator[,2] <= f$Comp.Indicator[,3]+1e-10))
# Independently enumerate all 24 permutations for m=4; use OLS residuals.
perms <- function(v) {
  if (length(v)==1L) return(matrix(v,nrow=1))
  do.call(rbind, lapply(v,function(j) cbind(j,perms(setdiff(v,j)))))
}
P <- perms(1:4)
q <- numeric(4)
for (i in seq_len(nrow(P))) {
  order <- P[i,]
  for (k in seq_along(order)) {
    j <- order[k]
    if(k==1) q[j] <- q[j]+1 else {
      fit <- lm.fit(cbind(1,X[,order[seq_len(k-1L)],drop=FALSE]),X[,j])
      q[j] <- q[j]+sum(fit$residuals^2)/sum((X[,j]-mean(X[,j]))^2)
    }
  }
}
stopifnot(near(q/nrow(P), f$q))
# Single indicator, polarity, and invalid-input checks.
stopifnot(near(dlp_rw(Z[,1,drop=FALSE],use_normat=FALSE)$Weights,1))
N <- dlp_rw_normalize(cbind(A=c(2,4),B=c(10,20)),polarity="A")
stopifnot(near(N$Z,c(0,1,1,0)))
stopifnot(inherits(try(dlp_rw(cbind(Z,C=1),use_normat=FALSE),silent=TRUE),"try-error"))
stopifnot(inherits(try(dlp_rw(Z,p=.5,use_normat=FALSE),silent=TRUE),"try-error"))
cat("All DLp-RW tests passed.\n")
