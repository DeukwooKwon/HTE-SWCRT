library(gee)
source("GeneratingData.R")

sig.level <- .05
dl <- DataGenerating(# for parameter settings, please refer to the function DataGenerating; for empirical type-I error, please make sure the last element of bete is 0
                     )
e <- dl$e
W <- dl$W
X <- dl$X
cluster.label <- dl$Cluster
sig.bin <- sapply(seq_len(ncol(dl$Y)), function(i){
  Y <- dl$Y[,i]
  fit <- gee(Y ~ factor(e) + W * X, id = cluster.label, family = "binomial"(link = "logit"), corstr = "exchangeable", maxiter = 50, tol = 1e-4)
  output <- summary(fit)
  return(1 - pnorm(abs(output$coefficients["W:X", 3])) < (sig.level / 2))
})
(emp.power <- mean(sig.bin))