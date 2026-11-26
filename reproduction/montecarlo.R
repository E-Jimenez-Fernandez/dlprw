# Run from this directory, or use Rscript montecarlo.R.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) setwd(dirname(normalizePath(sub('^--file=', '', file_arg[1]))))
library(dlprw)

B <- 100L  # Increase to 1000L for the final paper.
NS <- c(100L, 500L)
P <- c(1, 2, 3, 5, 10)
set.seed(20261126)
scenarios <- c('independent','moderate','strong','duplicate','correlated_duplicate','nonlinear')

# DP2: classical standard-deviation scaling or matched min-max input.
dp2 <- function(X, classic = TRUE, tol = 1e-8, maxit = 100L) {
  X <- as.matrix(X)
  D <- if (classic) sweep(sweep(X, 2, apply(X, 2, min), '-'), 2, apply(X, 2, sd), '/') else X
  V <- sweep(D, 2, colMeans(D), '-'); U <- sweep(V, 2, sqrt(colSums(V^2)), '/')
  m <- ncol(D); CI <- rowSums(D); q <- rep(1, m); seen <- character(); status <- 'maxit'
  residual <- function(j,S) {
    if (!length(S)) return(1)
    z <- svd(U[,S,drop=FALSE], nv=0)
    Q <- z$u[,z$d > 1e-10*z$d[1],drop=FALSE]
    r <- U[,j] - drop(Q %*% crossprod(Q,U[,j]))
    min(1,max(0,sum(r^2)))
  }
  ordering <- function(y) {
    rho <- abs(drop(cor(D,y)))
    if (anyNA(rho)) stop('Undefined DP2 correlation')
    order(-rho,seq_len(m))
  }
  for (iter in seq_len(maxit)) {
    ord <- ordering(CI); key <- paste(ord,collapse='-')
    if (key %in% seen) {status <- 'cycle'; break}
    seen <- c(seen,key); q <- rep(1,m)
    if (m>1) for (k in 2:m) q[ord[k]] <- residual(ord[k],ord[seq_len(k-1L)])
    new <- drop(D %*% q)
    if (abs(mean(new)-mean(CI)) <= tol) {CI <- new; status <- 'mean_tolerance'; break}
    CI <- new
    if (identical(ord,ordering(CI))) {status <- 'fixed_order'; break}
  }
  list(scores=CI,weights=q/sum(q),status=status)
}

kt <- function(x,y) cor(x,y,method='kendall')
fit_all <- function(X) {
  fit <- dlp_rw(X,p=P)
  list(fit=fit, C=fit$Comp.Indicator, d=dp2(X), dm=dp2(fit$normat,FALSE))
}
make_data <- function(n) {
  a <- runif(n); u <- runif(n); c <- runif(n); d <- runif(n)
  A <- cbind(X1=a,X2=u,X3=c,X4=d)
  B <- cbind(X1=a,X2=.5*a+.5*u,X3=c,X4=d)
  C <- cbind(X1=a,X2=.9*a+.1*u,X3=c,X4=d)
  setNames(list(A,B,C,cbind(A,X5=a),cbind(C,X5=a),
                cbind(X1=a,X2=(2*a-1)^2,X3=c,X4=d)),scenarios)
}
W <- list(); M <- list(); D <- list(); T <- list(); i <- 0L; k <- 0L
for (n in NS) for (b in seq_len(B)) {
  Xs <- make_data(n); keep <- sort(sample.int(n,floor(.9*n)))
  full <- lapply(Xs,fit_all)
  for (s in scenarios) {
    X <- Xs[[s]]; a <- full[[s]]; r <- fit_all(X[keep,,drop=FALSE]); f <- a$fit
    i <- i+1L
    W[[i]] <- data.frame(n=n,rep=b,scenario=s,indicator=seq_len(ncol(X)),q=unname(f$q),weight=unname(f$Weights),dp2_q_norm=a$d$weights,dp2_matched_weight=a$dm$weights)
    T[[i]] <- data.frame(n=n,rep=b,scenario=s,status_dp2=a$d$status,status_matched=a$dm$status,status_delete_dp2=r$d$status,status_delete_matched=r$dm$status)
    M[[i]] <- do.call(rbind,lapply(seq_along(P),function(j) data.frame(n=n,rep=b,scenario=s,p=P[j],tau_p1=kt(a$C[,1],a$C[,j]),tau_dp2=kt(a$C[,j],a$d$scores),tau_dp2_matched=kt(a$C[,j],a$dm$scores),tau_delete=kt(a$C[keep,j],r$C[,j]),tau_delete_dp2=kt(a$d$scores[keep],r$d$scores),tau_delete_matched=kt(a$dm$scores[keep],r$dm$scores),mean_score=mean(a$C[,j]),mean_gap=mean(a$C[,j]-a$C[,1]),mean_abs_weight_change=mean(abs(f$Weights-r$fit$Weights)),status_dp2=a$d$status,status_matched=a$dm$status)))
    if (s %in% c('duplicate','correlated_duplicate')) {
      base <- full[[if(s=='duplicate') 'independent' else 'strong']]; k <- k+1L
      gchange <- function(w,v) w[1]+w[length(w)]-v[1]
      D[[k]] <- do.call(rbind,lapply(seq_along(P),function(j) data.frame(n=n,rep=b,scenario=s,p=P[j],tau_dup=kt(base$C[,j],a$C[,j]),max_abs_score_change=max(abs(base$C[,j]-a$C[,j])),group_weight_change=gchange(f$Weights,base$fit$Weights),tau_dup_dp2=kt(base$d$scores,a$d$scores),tau_dup_matched=kt(base$dm$scores,a$dm$scores),group_change_dp2_q=gchange(a$d$weights,base$d$weights),group_change_matched=gchange(a$dm$weights,base$dm$weights))))
    }
  }
  if(b%%10L==0L) message('n=',n,', replication ',b,'/',B)
}
W <- do.call(rbind,W); M <- do.call(rbind,M); D <- do.call(rbind,D); T <- do.call(rbind,T)
summarize <- function(d,keys) {
  nums <- setdiff(names(d)[vapply(d,is.numeric,logical(1))],c(keys,'rep'))
  groups <- split(d,interaction(d[keys],drop=TRUE))
  do.call(rbind,lapply(groups,function(z) {
    out <- z[1,keys,drop=FALSE]
    for (v in nums) {out[[v]] <- mean(z[[v]]); out[[paste0(v,'_sd')]] <- sd(z[[v]]); out[[paste0(v,'_q05')]] <- unname(quantile(z[[v]],.05)); out[[paste0(v,'_q95')]] <- unname(quantile(z[[v]],.95))}
    out
  }))
}
dir.create('R_results',showWarnings=FALSE)
for (nm in c('W','M','D')) {
  z <- get(nm); name <- switch(nm,W='weights',M='metrics',D='duplication'); keys <- if(nm=='W') c('n','scenario','indicator') else c('n','scenario','p')
  write.csv(z,file.path('R_results',paste0(name,'_raw.csv')),row.names=FALSE)
  write.csv(summarize(z,keys),file.path('R_results',paste0(name,'_summary.csv')),row.names=FALSE)
}
write.csv(T,'R_results/dp2_status.csv',row.names=FALSE)
# Base R plots (R results, not the independent Python replication)
ms <- summarize(M,c('n','scenario','p')); z <- ms[ms$n==max(NS),]
pdf('R_results/kendall.pdf',width=9,height=6)
# Direct indexing avoids depending on ordering of aggregated rows.
mat <- sapply(scenarios,function(s) {v<-z[z$scenario==s,]; v$tau_p1[match(P,v$p)]})
matplot(P,mat,type='b',lty=1,pch=seq_along(scenarios),xlab='p',ylab='Mean Kendall tau vs p=1',ylim=c(.4,1))
legend('bottomleft',legend=scenarios,col=seq_along(scenarios),pch=seq_along(scenarios),lty=1,cex=.8)
dev.off()
print(table(T$status_dp2)); print(table(T$status_matched))
