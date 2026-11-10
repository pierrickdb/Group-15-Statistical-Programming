# =============================================================================
# Changepoint Estimation
#
# Authors:
# Pierrick de Boer: S3034691
# Fagan Ganbarli: S3015664
# Raneem Alluhaydan: S3044929
#
# Contributions:
# Fagan Ganbarli (%): Task 1-3
# Raneem Alluhaydan (%): Task 4-5
# Pierrick de Boer (%): Task 6-7
#
# Overview:
# 
#
# Main functions:
# 1. 
# =============================================================================


# Task 1: Generate the dataset.

set.seed(12345)

n <- 100

time <- 1:n

tau_true <- c(1, 30, 70, n+1) 

m_true <- length(tau_true)-2 

theta_true <- c(10, 2, 7) 

mu_true <- numeric(n) 

for(i in seq_along(theta_true)) {
  mu_true[tau_true[i]:(tau_true[i+1]-1)] <- theta_true[i]
}

y <- rpois(n, lambda = mu_true) ## simulate data y

plot(time, y, 
     pch = 16, 
     col = "grey")

lines(mu_true, 
      col="black", 
      lwd=2, 
      type ="s")

abline(v = tau_true[-c(1,length(tau_true))], 
       lty="dashed")



# Task 2(a) : 
poisson_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  n <- length(y)
  
  X <- outer(1:n, 1:n, ">=") * 1
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  neg_logik <- -sum(dpois(y, lambda = mu, log = TRUE))
  
  penalty <- lambda * sum(sqrt(beta[-1] ^ 2 + eps))
  
  return(neg_logik + penalty)
}


# Task 2(b) : 
poisson_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  n <- length(y)
  
  X <- outer(1:n, 1:n, ">=") * 1
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  grad <- as.vector(t(X) %*% (mu - y))
  
  grad[-1] <- grad[-1] + lambda * beta[-1] / sqrt(beta[-1] ^ 2 + eps)
  
  return(grad)
}


# Task 2(c) : 
poisson_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  n <- length(y)
  
  X <- outer(1:n, 1:n, ">=") * 1
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  H <- crossprod(X, sweep(X, 1, mu, "*"))
  
  penalty_hess <- c(0, lambda * eps / (beta[-1] ^ 2 + eps) ^ (3/2))
  
  diag(H) <- diag(H) + penalty_hess
  
  return(H)
}


# Task 2(b) :
finite_grad <- function(fn, beta, h = 1e-5) {
  
  t <- length(beta)
  g <- numeric(t)
  
  for (i in 1:t) {
    
    beta_plus <- beta_minus <- beta
    
    beta_plus[i] <- beta_plus[i] + h
    beta_minus[i] <- beta_minus[i] - h
    
    g[i] <- (fn(beta_plus) - fn(beta_minus)) / (2 * h)
  }
  
  return(g)
}


finite_hess <- function(fn, beta, h = 1e-5) {
  
  t <- length(beta)
  H <- matrix(0, t, t)
  
  for (i in 1:t) {
    
    beta_plus <- beta_minus <- beta
    
    beta_plus[i] <- beta_plus[i] + h
    beta_minus[i] <- beta_minus[i] - h
    
    H[, i] <- (fn(beta_plus) - fn(beta_minus)) / (2 * h)
  }
  
  return(H)
}

# Task 2(d): Checking
beta_test <- c(log(mean(y)), rep(0.05, n - 1))

grad_analytical <- poisson_lasso_grad(beta_test, y)

grad_numerical <- finite_grad(
  function(b) poisson_lasso_objective(b, y),
  beta_test
)

max(abs(grad_analytical - grad_numerical))


hess_analytical <- poisson_lasso_hess(beta_test, y)

hess_numerical <- finite_hess(
  function(b) poisson_lasso_grad(b, y),
  beta_test
)

max(abs(hess_analytical - hess_numerical))


# Task 3: 
beta_start <- c(log(mean(y)), rep(0, n - 1))

fit <- optim(
  par = beta_start,
  fn = poisson_lasso_objective,
  gr = poisson_lasso_grad,
  y = y,
  lambda = 1,
  method = "BFGS",
  control = list(maxit = 10000, reltol = 1e-9)
)

fit$convergence
fit$message

beta_hat <- fit$par
X <- outer(1:n, 1:n, ">=") * 1
mu_hat <- as.vector(exp(X %*% beta_hat))

plot(time, y,
     pch = 16,
     col = "grey",
     xlab = "Time",
     ylab = "Mean / Observations",
     main = "Poisson Changepoint Model")

lines(time, mu_true,
      col = "black",
      lwd = 2,
      type = "s")

lines(time, mu_hat,
      col = "blue",
      lwd = 2,
      type = "s")

legend("topright",
       legend = c("Observation data",
                  "True mean",
                  "Estimated mean"),
       col = c("grey", "black", "blue"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2))

















