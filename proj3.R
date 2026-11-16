# =============================================================================
# Changepoint Estimation
#
#
# GitHub repository: https://github.com/pierrickdb/Group-15-Statistical-Programming
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
# 1. poisson_lasso_objective
# 2. poisson_lasso_grad
# 3. poisson_lasso_hess
# 4. finite_grad
# 5. finite_hess
# =============================================================================

# Function definitions

# Task 2(a): Objective function

poisson_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6, X) {
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  neg_logik <- -sum(dpois(y, lambda = mu, log = TRUE))
  
  penalty <- lambda * sum(sqrt(beta[-1] ^ 2 + eps))
  
  return(neg_logik + penalty)
}



# Task 2(b): Gradient function

poisson_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6, X) {
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  grad <- as.vector(t(X) %*% (mu - y))
  
  grad[-1] <- grad[-1] + lambda * beta[-1] / sqrt(beta[-1] ^ 2 + eps)
  
  return(grad)
}


# Task 2(c): Hessian function

poisson_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6, X) {
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  H <- crossprod(X, sweep(X, 1, mu, "*"))
  
  penalty_hess <- c(0, lambda * eps / (beta[-1] ^ 2 + eps) ^ (3/2))
  
  diag(H) <- diag(H) + penalty_hess
  
  return(H)
}


# Task 2: Numerical gradient

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

# Task 2: Numerical Hessian

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


# =============================================================================
# Analysis and testing
# =============================================================================

# Task 1: Generate the data next set


set.seed(12345)

n <- 100

X <- outer(1:n, 1:n, ">=") * 1 



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

# Task 2: Checking the gradient and Hessian

beta_test <- c(log(mean(y)), rep(0.05, n - 1))

# Check gradient
grad_analytical <- poisson_lasso_grad(beta_test, y, X = X)

grad_numerical <- finite_grad(
  function(b) poisson_lasso_objective(b, y, X = X),
  beta_test
)

max(abs(grad_analytical - grad_numerical))


# Check Hessian

hess_analytical <- poisson_lasso_hess(beta_test, y, X = X)

hess_numerical <- finite_hess(
  function(b) poisson_lasso_grad(b, y, X = X),
  beta_test
)

max(abs(hess_analytical - hess_numerical))

# Task 3: Fit the model with lambda = 1

beta_start <- c(log(mean(y)), rep(0, n - 1))

fit <- optim(
  par = beta_start,
  fn = poisson_lasso_objective,
  gr = poisson_lasso_grad,
  y = y,
  lambda = 1,
  X = X,
  method = "BFGS",
  control = list(maxit = 10000, reltol = 1e-9)
)

fit$convergence
fit$message

beta_hat <- fit$par
mu_hat <- as.vector(exp(X %*% beta_hat))

# Plot
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


# Task 4: Choosing lambda using mBIC

lsp <- seq(-1, 5, by = 0.05)
mbic <- edf <- numeric(length(lsp))

b0 <- c(log(mean(y)), rep(0, n - 1))

for (i in 1:length(lsp)) {
  
  lambda <- exp(lsp[i])
  
  fit <- optim(b0, poisson_lasso_objective,
               gr = poisson_lasso_grad,
               y = y, lambda = lambda, X = X,
               method = "BFGS",
               control = list(maxit = 2000, reltol = 1e-9))
  
  if (fit$convergence != 0) {
    warning(paste("No convergence at", lambda))
    mbic[i] <- Inf
    next
  }
  
  b <- fit$par
  
  H1 <- poisson_lasso_hess(b, y, lambda, X = X)
  H0 <- poisson_lasso_hess(b, y, 0, X = X)
  
  edf[i] <- sum(diag(solve(H1, H0)))
  
  mu <- exp(cumsum(b))
  ll <- sum(dpois(y, mu, log = TRUE))
  
  mbic[i] <- 3 * (edf[i] - 1) * log(n) - 2 * ll
}

i.opt <- which.min(mbic)
lambda_best <- exp(lsp[i.opt])

lambda_best

plot(lsp, mbic, type = "l",
     xlab = "log(lambda)", ylab = "mBIC")

points(lsp[i.opt], mbic[i.opt], col = 2, pch = 19)


# Task 5: Fit the final model

fit <- optim(b0, poisson_lasso_objective,
             gr = poisson_lasso_grad,
             y = y, lambda = lambda_best, X = X,
             method = "BFGS",
             control = list(maxit = 2000, reltol = 1e-9))

fit$convergence

b <- fit$par
mu_hat <- exp(cumsum(b))

# Plot

plot(time, y, pch = 16, col = "grey")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat, col = "red", lwd = 2, type = "s")

# Number of change points

H1 <- poisson_lasso_hess(b, y, lambda_best, X = X)
H0 <- poisson_lasso_hess(b, y, 0, X = X)

edf_best <- sum(diag(solve(H1, H0)))
m_hat <- floor(edf_best - 1)

# Change point locations

if (m_hat > 0) {
  tau_hat <- sort(order(abs(b[-1]),
                        decreasing = TRUE)[1:m_hat] + 1)
} else {
  tau_hat <- numeric(0)
}

abline(v = tau_hat, col = "blue", lty = 2)

# Final results

lambda_best
m_hat
tau_hat


# Task 6: Parametric bootstrap

B <- 200
mu_boot <- matrix(NA, nrow = B, ncol = n)   # row b = bootstrap estimate mu^(b)
conv_boot <- numeric(B)

set.seed(1)

for (j in 1:B) {
  
  y_sim <- rpois(n, lambda = mu_hat)        # simulate from best fitted model
  
  b_start <- c(log(mean(y_sim)), rep(0, n - 1))
  
  fit_sim <- optim(b_start, poisson_lasso_objective,
                   gr = poisson_lasso_grad,
                   y = y_sim, lambda = lambda_best, X = X,
                   method = "BFGS",
                   control = list(maxit = 10000, reltol = 1e-9))
  
  conv_boot[j] <- fit_sim$convergence
  
  if (fit_sim$convergence != 0) {
    warning(paste("No convergence in bootstrap sample", j))
    next
  }
  
  mu_boot[j, ] <- exp(cumsum(fit_sim$par))
}

# Number of non-converged fits (these are excluded below)
sum(conv_boot != 0)

# Average bootstrap mean function
mu_bar <- colMeans(mu_boot, na.rm = TRUE)

# Plot: data, true mean, estimate, bootstrap average
plot(time, y, pch = 16, col = "grey",
     xlab = "Time", ylab = "y",
     main = "Bootstrap average vs estimated mean")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat, col = "red", lwd = 2, type = "s")
lines(time, mu_bar, col = "darkgreen", lwd = 2, type = "s")
legend("topright",
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bootstrap average"),
       col = c("grey", "black", "red", "darkgreen"),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), lwd = c(NA, 2, 2, 2))

# Bias corrected mean function: mu_hat - (mu_bar - mu_hat)
mu_bc <- 2 * mu_hat - mu_bar

# 95% bootstrap quantiles at each time point
q_L <- apply(mu_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
q_U <- apply(mu_boot, 2, quantile, probs = 0.975, na.rm = TRUE)

# Bias corrected 95% limits: reflect the quantiles about mu_hat
# (the bootstrap spread mu**(b) - mu_hat stands in for mu_hat - mu)
q_L_bc <- pmax(2 * mu_hat - q_U, 0)   #a Poisson mean cannot be negative
q_U_bc <- 2 * mu_hat - q_L


# Task 7: Final plot

ylim <- c(0, 1.35 * max(c(y, mu_true, mu_hat, mu_bc, q_U_bc)))  # room for legend

plot(time, y, pch = 16, col = "grey",
     ylim = ylim,
     xlab = "Time k", ylab = expression("Count " * y[k] * " / mean " * mu[k]),
     main = "Poisson changepoint model: bias corrected fit with 95% CI")

polygon(c(time, rev(time)), c(q_L_bc, rev(q_U_bc)),
        col = adjustcolor("blue", alpha.f = 0.15), border = NA)

lines(time, q_L_bc, col = "blue", lty = 2, lwd = 1.5, type = "s")
lines(time, q_U_bc, col = "blue", lty = 2, lwd = 1.5, type = "s")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat, col = "red", lwd = 2, type = "s")
lines(time, mu_bc, col = "blue", lwd = 2, type = "s")

legend("topright",
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bias corrected mean", "Bias corrected 95% CI"),
       col = c("grey", "black", "red", "blue", "blue"),
       pch = c(16, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 1, 2),
       lwd = c(NA, 2, 2, 2, 1.5),
       bg = "white")


=======
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
# 1. poisson_lasso_objective
# 2. poisson_lasso_grad
# 3. poisson_lasso_hess
# 4. finite_grad
# 5. finite_hess
# =============================================================================

# Function definitions

# Task 2(a): Objective function

poisson_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6, X) {
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  neg_logik <- -sum(dpois(y, lambda = mu, log = TRUE))
  
  penalty <- lambda * sum(sqrt(beta[-1] ^ 2 + eps))
  
  return(neg_logik + penalty)
}



# Task 2(b): Gradient function

poisson_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6, X) {
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  grad <- as.vector(t(X) %*% (mu - y))
  
  grad[-1] <- grad[-1] + lambda * beta[-1] / sqrt(beta[-1] ^ 2 + eps)
  
  return(grad)
}


# Task 2(c): Hessian function

poisson_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6, X) {
  
  eta <- as.vector(X %*% beta)
  
  mu <- exp(eta)
  
  H <- crossprod(X, sweep(X, 1, mu, "*"))
  
  penalty_hess <- c(0, lambda * eps / (beta[-1] ^ 2 + eps) ^ (3/2))
  
  diag(H) <- diag(H) + penalty_hess
  
  return(H)
}


# Task 2: Numerical gradient

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

# Task 2: Numerical Hessian

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


# =============================================================================
# Analysis and testing
# =============================================================================

# Task 1: Generate the data next set


set.seed(12345)

n <- 100

X <- outer(1:n, 1:n, ">=") * 1 



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

# Task 2: Checking the gradient and Hessian

beta_test <- c(log(mean(y)), rep(0.05, n - 1))

# Check gradient
grad_analytical <- poisson_lasso_grad(beta_test, y, X = X)

grad_numerical <- finite_grad(
  function(b) poisson_lasso_objective(b, y, X = X),
  beta_test
)

max(abs(grad_analytical - grad_numerical))


# Check Hessian

hess_analytical <- poisson_lasso_hess(beta_test, y, X = X)

hess_numerical <- finite_hess(
  function(b) poisson_lasso_grad(b, y, X = X),
  beta_test
)

max(abs(hess_analytical - hess_numerical))

# Task 3: Fit the model with lambda = 1

beta_start <- c(log(mean(y)), rep(0, n - 1))

fit <- optim(
  par = beta_start,
  fn = poisson_lasso_objective,
  gr = poisson_lasso_grad,
  y = y,
  lambda = 1,
  X = X,
  method = "BFGS",
  control = list(maxit = 10000, reltol = 1e-9)
)

fit$convergence
fit$message

beta_hat <- fit$par
mu_hat <- as.vector(exp(X %*% beta_hat))

# Plot
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


# Task 4: Choosing lambda using mBIC

lsp <- seq(-1, 5, by = 0.05)
mbic <- edf <- numeric(length(lsp))

b0 <- c(log(mean(y)), rep(0, n - 1))

for (i in 1:length(lsp)) {
  
  lambda <- exp(lsp[i])
  
  fit <- optim(b0, poisson_lasso_objective,
               gr = poisson_lasso_grad,
               y = y, lambda = lambda, X = X,
               method = "BFGS",
               control = list(maxit = 2000, reltol = 1e-9))
  
  if (fit$convergence != 0) {
    warning(paste("No convergence at", lambda))
    mbic[i] <- Inf
    next
  }
  
  b <- fit$par
  
  H1 <- poisson_lasso_hess(b, y, lambda, X = X)
  H0 <- poisson_lasso_hess(b, y, 0, X = X)
  
  edf[i] <- sum(diag(solve(H1, H0)))
  
  mu <- exp(cumsum(b))
  ll <- sum(dpois(y, mu, log = TRUE))
  
  mbic[i] <- 3 * (edf[i] - 1) * log(n) - 2 * ll
}

i.opt <- which.min(mbic)
lambda_best <- exp(lsp[i.opt])

lambda_best

plot(lsp, mbic, type = "l",
     xlab = "log(lambda)", ylab = "mBIC")

points(lsp[i.opt], mbic[i.opt], col = 2, pch = 19)


# Task 5: Fit the final model

fit <- optim(b0, poisson_lasso_objective,
             gr = poisson_lasso_grad,
             y = y, lambda = lambda_best, X = X,
             method = "BFGS",
             control = list(maxit = 2000, reltol = 1e-9))

fit$convergence

b <- fit$par
mu_hat <- exp(cumsum(b))

# Plot

plot(time, y, pch = 16, col = "grey")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat, col = "red", lwd = 2, type = "s")

# Number of change points

H1 <- poisson_lasso_hess(b, y, lambda_best, X = X)
H0 <- poisson_lasso_hess(b, y, 0, X = X)

edf_best <- sum(diag(solve(H1, H0)))
m_hat <- floor(edf_best - 1)

# Change point locations

if (m_hat > 0) {
  tau_hat <- sort(order(abs(b[-1]),
                        decreasing = TRUE)[1:m_hat] + 1)
} else {
  tau_hat <- numeric(0)
}

abline(v = tau_hat, col = "blue", lty = 2)

# Final results

lambda_best
m_hat
tau_hat


