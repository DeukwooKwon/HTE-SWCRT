require(boot)
require(MASS)

I.matrix <- function(n) diag(rep(1, n))
J.matrix <- function(n) matrix(1, nrow = n, ncol = n)

DataGenerating <- function(m = 50, # Cluster size
                           n = 15, # number of cases with covariate = 1 within each cluster (m * prevalence)
                           beta = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), 2), # fixed effects
                           alloc = c(1, 1, 2, 2, 3, 3, 4, 4), # Design info: time transitioning to intervention for each cluster
                           ICC = 0.1,
                           CAC = 0.8,
                           n.sims = 1000, # number of simulated samples
                           seed.value = 12345678 # set the seed
                           ){
  I <- length(alloc)
  cluster.index <- 1:I
  S <- length(unique(alloc))
  T.period <- S + 1
  if (T.period != length(beta) - 3){stop("Number of time fixed effect does not match the trial design!")}
  X.list <- lapply(cluster.index, function(i){
    rep(rep(c(1,0), c(n, (m - n))), T.period)
  })
  X <- unlist(X.list)
  state.ini <- rep(0, I)
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
  
  # Specifying the marginal means and correlation structure
  u.list <- lapply(M.list, function(i){
    inv.logit(as.numeric(i %*% beta))
  })
  Lsq.list <- lapply(u.list, function(i){
    diag(sqrt(i * (1 - i)))
  })
  # labeling correlated individuals for the GLMM
  cluster.label <- rep(cluster.index, each = (m * T.period))
  # period.label <- rep(rep(1:T.period, each = m), 8)
  u <- inv.logit(M %*% beta)
  R.list <- lapply(rep(m, I), function(i){
    ICC * CAC * J.matrix(i * T.period) + ICC * (1 - CAC) * kronecker(I.matrix(T.period), J.matrix(i)) + (1 - ICC) * I.matrix(i * T.period)
  })
  A.list <- lapply(cluster.index, function(i){
    Lsq.list[[i]] %*% R.list[[i]] %*% Lsq.list[[i]] + u.list[[i]] %*% t(u.list[[i]])
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
  
  ### Correction in case some of the calculated correlation matrices are not positively definite
  Rho.PD.list <- lapply(cluster.index, function(i){
    Eig.original <- eigen(Rho.list[[i]])
    Eig.new <- ifelse(Eig.original$values < 0, 0, Eig.original$values)
    newMat <- Eig.original$vectors %*% diag(Eig.new) %*% t(Eig.original$vectors)
    newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
    return(newMat)
  })
  
  set.seed(seed.value)
  Y.mat <- replicate(n.sims, {
    Y.list <- lapply(cluster.index, function(i){
      mvrnorm(n = 1, mu = rep(0, length(u.list[[i]])), Sigma = Rho.PD.list[[i]]) < qnorm(u.list[[i]])
    })
    Y <- unlist(Y.list)
  })
  return(list("e" = e, "W" = W, "X" = X, "Cluster" = cluster.label, "Y" = Y.mat))
}