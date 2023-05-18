library(VGAM)
library(boot)
library(MASS)
library(gee)
m <- 160
p <- 0.3
n <- 48 # n = m * p
beta <- c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.35), log(1.5), log(1.5))
beta0 <- c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.35), log(1.5), 0)
# Specifying the design matrix
alloc <- c(1, 1, 2, 2, 3, 3, 4, 4)
I <- length(alloc)
cluster.index <- 1:I
S <- length(unique(alloc))
T.period <- S + 1
X.list <- lapply(cluster.index, function(i){
  rep(rep(c(1,0), c(n, (m - n))), T.period)
})
X <- unlist(X.list)
state.ini <- rep(0, 8)
SW.design <- sapply(1:T.period, function(s){
  state.ini[which(alloc < s)] <- 1
  return(state.ini)
})
# Specifying the regression model matrix
W.list <- lapply(cluster.index, function(i){
  rep(SW.design[i,], each = m)
})
W <- unlist(W.list)
e.list <- lapply(cluster.index, function(i){
  rep(0:S, each = m)
})
e <- unlist(e.list)
M.list <- lapply(cluster.index, function(i){
  model.matrix(~ factor(e.list[[i]]) + W.list[[i]] * X.list[[i]])
})
M <- model.matrix(~ factor(e) + W * X)
u.list <- lapply(M.list, function(i){
  inv.logit(as.numeric(i %*% beta))
})
u0.list <- lapply(M.list, function(i){
  inv.logit(as.numeric(i %*% beta0))
})
Lsq.list <- lapply(u.list, function(i){
  diag(sqrt(i * (1 - i)))
})
Lsq0.list <- lapply(u0.list, function(i){
  diag(sqrt(i * (1 - i)))
})
# labeling correlated individuals for the GLMM
cluster.label <- rep(cluster.index, each = (m * T.period))
# period.label <- rep(rep(1:T.period, each = m), 8)
# Significance level
sig.level <- 0.05
u <- inv.logit(M %*% beta)
u0 <- inv.logit(M %*% beta0)
ICC <- 0.1
I.matrix <- function(n) diag(rep(1, n))
J.matrix <- function(n) matrix(1, nrow = n, ncol = n)
R.list <- lapply(rep(m, I), function(i){
  ICC * J.matrix(i * T.period) + (1 - ICC) * I.matrix(i * T.period)
})
A.list <- lapply(cluster.index, function(i){
  Lsq.list[[i]] %*% R.list[[i]] %*% Lsq.list[[i]] + u.list[[i]] %*% t(u.list[[i]])
})
A0.list <- lapply(cluster.index, function(i){
  Lsq0.list[[i]] %*% R.list[[i]] %*% Lsq0.list[[i]] + u0.list[[i]] %*% t(u0.list[[i]])
})
Rho.list <- lapply(cluster.index, function(i){
  outer(X = 1:length(u.list[[i]]), Y = 1:length(u.list[[i]]), FUN = Vectorize(function(j,k){
    if (k == j) {return(1)} else{
      x <- u.list[[i]][j]
      y <- u.list[[i]][k]
      if (R.list[[i]][j,k] == 0) {return(0)} else{
        a <- A.list[[i]][j,k]
        f <- function(r) pbinormcop(x, y, rho = r) - a
        rho <- uniroot(f, interval = c(-0.99999,0.99999))
        return(rho$root)
      }
    }
  }))
})
Rho0.list <- lapply(cluster.index, function(i){
  outer(X = 1:length(u0.list[[i]]), Y = 1:length(u0.list[[i]]), FUN = Vectorize(function(j,k){
    if (k == j) {return(1)} else{
      x <- u0.list[[i]][j]
      y <- u0.list[[i]][k]
      if (R.list[[i]][j,k] == 0) {return(0)} else{
        a <- A0.list[[i]][j,k]
        f <- function(r) pbinormcop(x, y, rho = r) - a
        rho <- uniroot(f, interval = c(-0.99999,0.99999))
        return(rho$root)
      }
    }
  }))
})
Rho.PD.list <- lapply(cluster.index, function(i){
  Eig.original <- eigen(Rho.list[[i]])
  Eig.new <- ifelse(Eig.original$values < 0, 0, Eig.original$values)
  newMat <- Eig.original$vectors %*% diag(Eig.new) %*% t(Eig.original$vectors)
  newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
  return(newMat)
})
Rho0.PD.list <- lapply(cluster.index, function(i){
  Eig.original <- eigen(Rho0.list[[i]])
  Eig.new <- ifelse(Eig.original$values < 0, 0, Eig.original$values)
  newMat <- Eig.original$vectors %*% diag(Eig.new) %*% t(Eig.original$vectors)
  newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
  return(newMat)
})
n.sims <- 1000

set.seed(12345678)
simu.sig <- replicate(n.sims, {
  Y.list <- lapply(cluster.index, function(i){
    mvrnorm(n = 1, mu = rep(0, length(u.list[[i]])), Sigma = Rho.PD.list[[i]]) < qnorm(u.list[[i]])
  })
  Y <- unlist(Y.list)
  fit <- gee(Y ~ factor(e) + W * X, id = cluster.label, family = "binomial"(link = "logit"), corstr = "exchangeable", maxiter = 50, tol = 1e-4)
  output <- summary(fit)
  Y0.list <- lapply(cluster.index, function(i){
    mvrnorm(n = 1, mu = rep(0, length(u0.list[[i]])), Sigma = Rho0.PD.list[[i]]) < qnorm(u0.list[[i]])
  })
  Y0 <- unlist(Y0.list)
  fit0 <- gee(Y0 ~ factor(e) + W * X, id = cluster.label, family = "binomial"(link = "logit"), corstr = "exchangeable", maxiter = 50, tol = 1e-4)
  output0 <- summary(fit0)
  c(1 - pnorm(abs(output$coefficients["W:X", 3])) < (sig.level / 2),
    1 - pnorm(abs(output0$coefficients["W:X", 3])) < (sig.level / 2))
})
rowMeans(simu.sig)
save(simu.sig, file = "~/n3out/GEE_MD3_Medium.RData")
