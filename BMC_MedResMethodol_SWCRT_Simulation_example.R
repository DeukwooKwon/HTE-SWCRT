library(boot)
library(lme4)
library(optimx)

### Specific SW-CRT design
# Input: vector of numbers of individuals within each cluster per period
m <- rep(120, 8)
cluster.index <- 1:length(m)

# Intervention transitioning allocation (represented by indices of clusters, 2 transitions per period)
alloc <- c(1, 1, 2, 2, 3, 3, 4, 4)
S <- length(unique(alloc))
T.period <- S + 1
state.ini <- rep(0, length(m))
SW.design <- sapply(1:T.period, function(s){
  state.ini[which(alloc < s)] <- 1
  return(state.ini)
})

# Model matrices
W.list <- lapply(cluster.index, function(i){
  rep(SW.design[i,], each = m[i])
})
W <- unlist(W.list)
e.list <- lapply(cluster.index, function(i){
  rep(0:S, each = m[i])
})
e <- unlist(e.list)
cluster.label <- rep(cluster.index, (m * T.period))
period.label <- unlist(lapply(m, function(i){
  rep(1:T.period, each = i)
}))
individual.label <- unlist(lapply(m, function(i){
  rep(1:i, T.period)
}))

# Model (true) parameters
# linear coefficients for scenarios (varying interaction coefficient)

# beta <- c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), 0)
 beta <- c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(1.5))
# beta <- c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2))


# significant level
sig.level <- 0.05

# standard deviations of random effects
sigma.alpha <- 1
sigma.nu <- 0.5
sigma.psi <- 0.75

### Number of simulated samples
n.sims <- 1000

# Prevalence level scenarios (Beta-binomial data generating process for binary covariate)
set.seed(12345678)
# X.cluster.means <- rbeta(length(m), shape1 = 0.9, shape2 = 2.1) # Expected prevalence: 0.3
X.cluster.means <- rbeta(length(m), shape1 = 2.1, shape2 = 2.1) # Expected prevalence: 0.5
# X.cluster.means <- rbeta(length(m), shape1 = 4.9, shape2 = 2.1) # Expected prevalence: 0.7
X.list <- lapply(cluster.index, function(i){
  rep(rbinom(m[i], size = 1, prob = X.cluster.means[i]), T.period)
})
X <- unlist(X.list)
M.list <- lapply(cluster.index, function(i){
  model.matrix(~ factor(e.list[[i]]) + W.list[[i]] * X.list[[i]])
})
M <- model.matrix(~ factor(e) + W * X)

### Simulation (each data is generated using the codes within the function replicate)

simu.sig <- replicate(n.sims, {
  alpha <- rep(rnorm(length(m), mean = 0, sd = sigma.alpha), (m * T.period))
  nu <- rep(rnorm(T.period * length(m), mean = 0, sd = sigma.nu), rep(m, each = T.period))
  psi <- unlist(lapply(m, function(i){
    rep(rnorm(i, mean = 0, sd = sigma.psi), T.period)
  }))
  linear.predicator <- M %*% beta + alpha + nu + psi
  u <- inv.logit(linear.predicator)
  Y <- rbinom(length(u), size = 1, prob = u)
  fit <- glmer(Y ~ factor(e) + W * X + (1 | cluster.label) + (1 | cluster.label : period.label) + (1 | cluster.label : individual.label),
               family = "binomial",
               control = glmerControl(optimizer = "Nelder_Mead",
                                      boundary.tol = 0,
                                      check.scaleX = "ignore",
                                      check.conv.grad = .makeCC(action = "ignore", tol = 1e-4),
                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
  output <- summary(fit)
  return(output$coefficients["W:X", 4] < sig.level)
  
})
mean(simu.sig)
