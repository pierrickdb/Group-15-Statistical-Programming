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
# Fagan Ganbarli (33.33%): Task 1-3
# Raneem Alluhaydan (33.33%): Task 4-5
# Pierrick de Boer (33.33%): Task 6-7
#
# Overview:
# Estimate changepoints in Poisson data using a smoothed LASSO penalty.
# Choose lambda by minimising mBIC.
# Use parametric bootstrapping for bias correction and confidence intervals.
#
# Main functions:
# 1. poisson_lasso_objective: Calculates the penalised cost.
# 2. poisson_lasso_grad: Calculates the gradient.
# 3. poisson_lasso_hess: Calculates the Hessian.
# 4. fd_jacobian: Checks derivatives numerically.
# 5. fit_poisson_lasso: Fits the model using BFGS.
# 6. poisson_lasso_edf: Calculates effective degrees of freedom.
# =============================================================================

# Function definitions

# Task 2(a): Objective function

poisson_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  mu <- exp(cumsum(beta))
  
  neg_loglik <- -sum(dpois(y, lambda = mu, log = TRUE))
  
  penalty <- lambda * sum(sqrt(beta[-1]^2 + eps))
  
  return(neg_loglik + penalty)
}


# Task 2(b): Gradient function

poisson_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  mu <- exp(cumsum(beta))
  
  grad <- rev(cumsum(rev(mu - y)))
  
  grad[-1] <- grad[-1] +
    lambda * beta[-1] / sqrt(beta[-1]^2 + eps)
  
  return(grad)
}


# Task 2(c): Hessian function

poisson_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  n <- length(beta)
  mu <- exp(cumsum(beta))
  
  s <- rev(cumsum(rev(mu)))
  
  idx <- outer(1:n, 1:n, pmax)
  H <- matrix(s[idx], n, n)
  
  penalty_hess <- c(0, lambda * eps /
                      (beta[-1]^2 + eps)^(3/2))
  
  diag(H) <- diag(H) + penalty_hess
  
  return(H)
}


# Numerical differentiation
# Applied to the objective to check the gradient,
# and to the gradient to check the Hessian.

fd_jacobian <- function(fn, beta, h = 1e-5) {
  
  n <- length(beta)
  J <- matrix(0, length(fn(beta)), n)
  
  for (i in 1:n) {
    beta_plus <- beta_minus <- beta
    beta_plus[i] <- beta_plus[i] + h
    beta_minus[i] <- beta_minus[i] - h
    
    J[, i] <- (fn(beta_plus) - fn(beta_minus)) / (2 * h)
  }
  
  return(J)
}


# Fit the Poisson LASSO model

fit_poisson_lasso <- function(y, lambda = 1, beta_start = NULL) {
  
  if (is.null(beta_start)) {
    beta_start <- c(log(mean(y)), rep(0, length(y) - 1))
  }
  
  fit <- optim(
    beta_start, poisson_lasso_objective,
    gr = poisson_lasso_grad,
    y = y, lambda = lambda,
    method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-9)
  )
  
  if (fit$convergence != 0) {
    warning(paste("No convergence at lambda =", lambda))
  }
  
  return(fit)
}


# Calculate effective degrees of freedom

poisson_lasso_edf <- function(beta, y, lambda) {
  
  H1 <- poisson_lasso_hess(beta, y, lambda)
  H0 <- poisson_lasso_hess(beta, y, 0)
  
  edf <- sum(diag(solve(H1, H0)))
  
  return(edf)
}



# =============================================================================
# Analysis and testing
# =============================================================================

# Task 1: Generate the data set

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
     col = "grey",
     xlab = "Time",
     ylab = "Count",
     main = "Simulated Poisson Changepoint Data")

lines(mu_true, 
      col="black", 
      lwd=2, 
      type ="s")

abline(v = tau_true[-c(1,length(tau_true))], 
       lty="dashed")



# Task 2: Checking the gradient and Hessian

beta_test <- c(log(mean(y)), rep(0.05, n - 1))

grad_analytical <- poisson_lasso_grad(beta_test, y)

grad_numerical <- as.vector(fd_jacobian(
  function(b) poisson_lasso_objective(b, y),
  beta_test
))

hess_analytical <- poisson_lasso_hess(beta_test, y)

hess_numerical <- fd_jacobian(
  function(b) poisson_lasso_grad(b, y),
  beta_test
)

max(abs(grad_analytical - grad_numerical))
max(abs(hess_analytical - hess_numerical))

# Gradient error = 4.230387e-07
# Hessian error = 5.817456e-07
# Both are close to zero, confirming the analytical derivatives.


# Task 3: Fit the model with lambda = 1

beta_start <- c(log(mean(y)), rep(0, n - 1))

fit <- fit_poisson_lasso(y, lambda = 1, beta_start = beta_start)


fit$convergence
fit$message

# The model converged with lambda = 1.
# The fitted mean follows the trend but also fits some noise.

beta_hat <- fit$par
mu_hat1 <- exp(cumsum(beta_hat))

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

lines(time, mu_hat1, col = "blue", lwd = 2, type = "s")

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
mbic <- rep(Inf, length(lsp))
edf <- rep(NA_real_, length(lsp))

b0 <- beta_start

for (i in seq_along(lsp)) {
  
  lambda <- exp(lsp[i])
  
  fit_lambda <- fit_poisson_lasso(y, lambda, b0)
  
  if (fit_lambda$convergence != 0) {
    next
  }
  
  beta_lambda <- fit_lambda$par
  
  edf[i] <- poisson_lasso_edf(beta_lambda, y, lambda)
  
  mu_lambda <- exp(cumsum(beta_lambda))
  ll <- sum(dpois(y, mu_lambda, log = TRUE))
  
  mbic[i] <- 3 * (edf[i] - 1) * log(n) - 2 * ll
  
  # Use previous estimates as starting values
  b0 <- beta_lambda
}

i_opt <- which.min(mbic)
lambda_best <- exp(lsp[i_opt])

lambda_best

# mBIC is minimised at log(lambda) = 3, giving lambda = 20.08554.

plot(lsp, mbic, type = "l",
     xlab = "log(lambda)",
     ylab = "mBIC",
     main = "Choosing lambda using mBIC")

points(lsp[i_opt], mbic[i_opt],
       col = "red", pch = 19)




# Task 5: Fit the final model

fit_best <- fit_poisson_lasso(
  y, lambda_best,
  c(log(mean(y)), rep(0, n - 1))
)

# Check convergence
fit_best$convergence

beta_best <- fit_best$par
mu_hat_best <- exp(cumsum(beta_best))

# Number of changepoints
edf_best <- poisson_lasso_edf(beta_best, y, lambda_best)
m_hat <- floor(edf_best - 1)

# Changepoint locations
if (m_hat > 0) {
  tau_hat <- sort(order(abs(beta_best[-1]),
                        decreasing = TRUE)[1:m_hat] + 1)
} else {
  tau_hat <- numeric(0)
}

# Plot
plot(time, y, pch = 16, col = "grey",
     xlab = "Time",
     ylab = "Count / Mean",
     main = "Final Poisson Changepoint Model")

lines(time, mu_true, col = "black",
      lwd = 2, type = "s")

lines(time, mu_hat_best, col = "red",
      lwd = 2, type = "s")

abline(v = tau_hat, col = "blue", lty = 2)

legend("topright", cex = 0.55,
         legend = c("Data", "True mean",
                  "Estimated mean", "Changepoints"),
       col = c("grey", "black", "red", "blue"),
       pch = c(16, NA, NA, NA),
       lty = c(NA, 1, 1, 2),
       bty = "n")

# Final results
lambda_best
m_hat
tau_hat

# EDF = 5.864, giving 4 estimated changepoints.
# Estimated locations: 30, 69, 70 and 72.
# True locations: 30 and 70.
# The change near 70 is split into several estimated changepoints.
# The penalty shrinks the estimated mean changes.


# Task 6: Parametric bootstrap

# Bootstrap variation estimates uncertainty in the fitted mean.
# Reflecting the quantiles gives bias-corrected confidence limits.

B <- 200
mu_boot <- matrix(NA, nrow = B, ncol = n)  

conv_boot <- numeric(B)

set.seed(1)

for (j in 1:B) {
  
  y_sim <- rpois(n, lambda = mu_hat_best)       
  
  # simulate from best fitted model
  
  b_start <- c(log(mean(y_sim)), rep(0, n - 1))
  
  
  fit_sim <- fit_poisson_lasso(
    y_sim, lambda_best,
    beta_start = b_start
  )
  
  conv_boot[j] <- fit_sim$convergence
  
  #sanity check if not converges
  if (fit_sim$convergence != 0) {
    warning(paste("No convergence in bootstrap sample", j))
    next
  }
  
  mu_boot[j, ] <- exp(cumsum(fit_sim$par))
}

# Number of non-converged fits 
sum(conv_boot != 0)

# Average bootstrap mean function
mu_bar <- colMeans(mu_boot, na.rm = TRUE)

# Iterpretation of the average boot strap: 
# The bootstraped average may be further from the true mean
# because the simulations start from an already biased estimate.

# Plot: data, true mean, estimate, bootstrap average

plot(time, y, pch = 16, col = "grey",
     ylim = c(0, 1.3 * max(y, mu_true, mu_hat_best, mu_bar)),
     xlab = "Time", ylab = "Count / Mean",
     main = "Bootstrap average vs estimated mean")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat_best, col = "red", lwd = 2, type = "s")
lines(time, mu_bar, col = "darkgreen", lwd = 2, type = "s")
legend("topright",
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bootstrap average"),
       col = c("grey", "black", "red", "darkgreen"),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), lwd = c(NA, 2, 2, 2))

# Bias corrected mean function: mu_hat_best - (mu_bar - mu_hat_best)
mu_bc <- 2 * mu_hat_best - mu_bar

# 95% bootstrap quantiles at each time point
q_L <- apply(mu_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
q_U <- apply(mu_boot, 2, quantile, probs = 0.975, na.rm = TRUE)

# estimate the variability of the fitted mean.
q_L_bc <- pmax(2 * mu_hat_best - q_U, 0)
q_U_bc <- 2 * mu_hat_best - q_L


# Task 7: Final plot

ylim <- c(0, 22) 

plot(time, y, pch = 16, col = "grey",
     ylim = ylim,
     xlab = "Time k", ylab = expression("Count " * y[k] * " / mean " * mu[k]),
     main = "Bias-corrected Poisson model")

polygon(c(time, rev(time)), c(q_L_bc, rev(q_U_bc)),
        col = adjustcolor("blue", alpha.f = 0.15), border = NA)

lines(time, q_L_bc, col = "blue", lty = 2, lwd = 1.5, type = "s")
lines(time, q_U_bc, col = "blue", lty = 2, lwd = 1.5, type = "s")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_hat_best, col = "red", lwd = 2, type = "s")
lines(time, mu_bc, col = "blue", lwd = 2, type = "s")

legend("topright", cex = 0.7,
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bias corrected mean", "Bias corrected 95% CI"),
       col = c("grey", "black", "red", "blue", "blue"),
       pch = c(16, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 1, 2),
       lwd = c(NA, 2, 2, 2, 1.5),
       bg = "white")

# Interpretation of the final plot:
# Shrinkage pulls the segment means of mu_hat towards each other; the bias
# correction largely reverses this, especially in the middle segment.
round(c(mu_hat = mean(mu_hat_best[31:68]), mu_bc = mean(mu_bc[31:68])), 1)  

# The band is narrow within segments but widens sharply around k = 30 and k = 70:
# the bootstrap fits disagree on where the jump is, so mu_bar is smeared and
# mu_bc = 2*mu_hat - mu_bar overshoots just before the jump.
round(mu_bc[29], 1)

# Coverage: proportion of time points where the band contains the true mean
mean(q_L_bc <= mu_true & mu_true <= q_U_bc)
which(!(q_L_bc <= mu_true & mu_true <= q_U_bc))
# 0.97 is close to the nominal 95%. The misses are at the jump near k = 70,
# which the fit spreads over k = 69-72.
