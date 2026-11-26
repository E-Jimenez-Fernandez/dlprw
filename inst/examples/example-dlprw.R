library(dlprw)
set.seed(123)
n <- 100
income <- runif(n)
X <- data.frame(Income = income,
                Education = .7 * income + .3 * runif(n),
                Health = runif(n),
                Employment = runif(n))
result <- dlp_rw(X, polarity=1:4, p=c(1,2,3,5,10))
print(data.frame(Indicator=names(result$Weights), q=result$q,
                 Weight=result$Weights))
print(head(result$Comp.Indicator))
print(cor(result$Comp.Indicator,method="kendall"))
print(head(result$Ranking))
