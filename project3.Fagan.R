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

# FUNCTION DEFINITIONS

# Task 2(a): Objective function
# Calculate the penalised negative log-likelihood
# for the Poisson changepoint model.

poisson_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  # Inputs:
  #   beta: Parameter vector of length n.
  #   y: Observed data vector.
  #   lambda: LASSO penalty parameter. Default value is 1.
  #   eps: Smoothing parameter. Default value is 1e-6.
  
  # Output:
  #   The penalised negative log-likelihood.
  
  # Method:
  #   Calculate the Poisson mean function, evaluate the
  #   negative log-likelihood and add the smoothed LASSO penalty.
  
  # Calculate the Poisson mean function.
  mu <- exp(cumsum(beta))
  
  # Calculate the negative Poisson log-likelihood.
  neg_loglik <- -sum(dpois(y, lambda = mu, log = TRUE))
  
  # Apply the smoothed LASSO penalty, excluding the first parameter.
  penalty <- lambda * sum(sqrt(beta[-1]^2 + eps))
  
  # Return the penalised objective function.
  return(neg_loglik + penalty)
}


# Task 2(b): Gradient function
# Calculate the gradient of the penalised objective function.

poisson_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  # Inputs:
  #   beta: Parameter vector of length n.
  #   y: Observed data vector.
  #   lambda: LASSO penalty parameter. Default value is 1.
  #   eps: Smoothing parameter. Default value is 1e-6.
  
  # Output:
  #   Gradient vector containing the first partial derivatives.
  
  # Method:
  #   Calculate the likelihood gradient using reverse cumulative
  #   sums and add the gradient of the smoothed penalty.
  
  # Calculate the Poisson mean function.
  mu <- exp(cumsum(beta))
  
  # Calculate the gradient of the negative log-likelihood.
  grad <- rev(cumsum(rev(mu - y)))
  
  # Add the penalty gradient, excluding the first parameter.
  grad[-1] <- grad[-1] +
    lambda * beta[-1] / sqrt(beta[-1]^2 + eps)
  
  # Return the gradient vector.
  return(grad)
}


# Task 2(c): Hessian function
# Calculate the Hessian matrix of the penalised objective function.

poisson_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6) {
  
  # Inputs:
  #   beta: Parameter vector of length n.
  #   y: Observed data vector.
  #   lambda: LASSO penalty parameter. Default value is 1.
  #   eps: Smoothing parameter. Default value is 1e-6.
  
  # Output:
  #   An n by n Hessian matrix containing the second derivatives.
  
  # Method:
  #   Calculate the likelihood Hessian using reverse cumulative sums
  #   and add the second derivatives of the penalty.
  
  # Store the number of parameters.
  n <- length(beta)
  mu <- exp(cumsum(beta))
  
  # Calculate reverse cumulative sums of the Poisson means.
  s <- rev(cumsum(rev(mu)))
  
  # Find the maximum index for each pair of parameters.
  idx <- outer(1:n, 1:n, pmax)
  H <- matrix(s[idx], n, n)
  
  # Calculate the second derivatives of the penalty.
  penalty_hess <- c(0, lambda * eps /
                      (beta[-1]^2 + eps)^(3/2))
  
  # Add the penalty derivatives to the diagonal.
  diag(H) <- diag(H) + penalty_hess
  
  # Return the Hessian matrix.
  return(H)
}


# Task 2(d): Numerical differentiation
# Use central finite differences to check the analytical
# gradient and Hessian functions.

fd_jacobian <- function(fn, beta, h = 1e-5) {
  
  # Inputs:
  #   fn: Function to differentiate.
  #   beta: Parameter vector at which derivatives are evaluated.
  #   h: Finite difference step size. Default value is 1e-5.
  
  # Output:
  #   Numerical Jacobian matrix.
  
  # Store the number of parameters.
  n <- length(beta)
  J <- matrix(0, length(fn(beta)), n)
  
  # Calculate numerical derivatives for each parameter.
  for (i in 1:n) {
    beta_plus <- beta_minus <- beta
    beta_plus[i] <- beta_plus[i] + h
    beta_minus[i] <- beta_minus[i] - h
    
    # Apply the central finite difference approximation.
    J[, i] <- (fn(beta_plus) - fn(beta_minus)) / (2 * h)
  }
  
  # Return the numerical Jacobian.
  return(J)
}


# Task 3: Fit the Poisson LASSO model using optim.
# Estimate the beta parameters by minimising the penalised
# negative log-likelihood using the BFGS method.

fit_poisson_lasso <- function(y, lambda = 1, beta_start = NULL) {
  
  # Inputs:
  #   y: Observed data vector.
  #   lambda: LASSO penalty parameter. Default value is 1.
  #   beta_start: Optional initial parameter vector.
  
  # Output:
  #   Optimisation results containing estimated parameters
  #   and convergence information.
  
  # Method:
  #   Initialise the parameters and minimise the objective
  #   function using BFGS with the analytical gradient.
  
  # Set the initial parameter values if none are provided.
  if (is.null(beta_start)) {
    beta_start <- c(log(mean(y)), rep(0, length(y) - 1))
  }
  
  # Fit the model using the BFGS optimisation method.
  fit <- optim(
    beta_start, poisson_lasso_objective,
    gr = poisson_lasso_grad,
    y = y, lambda = lambda,
    method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-9)
  )
  
  # Check whether the optimisation algorithm converged.
  if (fit$convergence != 0) {
    warning(paste("No convergence at lambda =", lambda))
  }
  
  # Return the optimisation results.
  return(fit)
}


# Task 4: Calculate the effective degrees of freedom.

# Calculate the effective degrees of freedom using
# the penalised and unpenalised Hessian matrices.
poisson_lasso_edf <- function(beta, y, lambda) {
  
  # Inputs:
  #   beta: Estimated parameter vector.
  #   y: Observed data vector.
  #   lambda: LASSO penalty parameter.
  
  # Output:
  #   Effective degrees of freedom of the fitted model.
  
  # Method:
  #   Evaluate the penalised and unpenalised Hessians
  #   and calculate the trace of their matrix product.
  
  # Calculate the penalised Hessian.
  H1 <- poisson_lasso_hess(beta, y, lambda)
  
  # Calculate the unpenalised Hessian.
  H0 <- poisson_lasso_hess(beta, y, 0)
  
  # Calculate the effective degrees of freedom.
  edf <- sum(diag(solve(H1, H0)))
  
  # Return the effective degrees of freedom.
  return(edf)
}

# =============================================================================
# Analysis and testing
# =============================================================================

# Task 1: Generate the data set

# Define the number of observations and time points.
set.seed(12345)
n <- 100
time <- 1:n

# Specify the true changepoint locations.
tau_true <- c(1, 30, 70, n + 1)

# Calculate the true number of changepoints.
m_true <- length(tau_true) - 2

# Define the Poisson mean for each segment.
theta_true <- c(10, 2, 7)

mu_true <- numeric(n)

# Assign the corresponding mean to each time segment.
for (i in seq_along(theta_true)) {
  mu_true[tau_true[i]:(tau_true[i + 1] - 1)] <- theta_true[i]
}

# Generate observations from the Poisson changepoint model.
y <- rpois(n, lambda = mu_true)

# Plot the simulated observations and true mean function.
plot(time, y,
     pch = 16,
     col = "grey",
     xlab = "Time",
     ylab = "Count",
     main = "Simulated Poisson Changepoint Data")

lines(time, mu_true,
      col = "black",
      lwd = 2,
      type = "s")

abline(v = tau_true[-c(1, length(tau_true))],
       lty = "dashed")


# Task 2: Checking the gradient and Hessian

# Choose test parameter values for numerical differentiation.
beta_test <- c(log(mean(y)), rep(0.05, n - 1))

# Calculate the analytical gradient.
grad_analytical <- poisson_lasso_grad(beta_test, y)

# Approximate the gradient using central finite differences.
grad_numerical <- as.vector(fd_jacobian(
  function(b) poisson_lasso_objective(b, y),
  beta_test
))

# Calculate the analytical Hessian.
hess_analytical <- poisson_lasso_hess(beta_test, y)

# Approximate the Hessian using finite differences of the gradient.
hess_numerical <- fd_jacobian(
  function(b) poisson_lasso_grad(b, y),
  beta_test
)

# Compare the analytical and numerical derivatives.
max(abs(grad_analytical - grad_numerical))
max(abs(hess_analytical - hess_numerical))

# Gradient error = 4.230387e-07
# Hessian error = 5.817456e-07
# Both are close to zero, confirming the analytical derivatives.


# Task 3: Fit the model with lambda = 1

# Initialise the parameters using the sample mean.
beta_start <- c(log(mean(y)), rep(0, n - 1))

# Fit the Poisson LASSO model using BFGS.
fit <- fit_poisson_lasso(y, lambda = 1, beta_start = beta_start)

# Check whether the optimisation algorithm converged.
fit$convergence
fit$message

# Stop if the optimisation did not converge.
if (fit$convergence != 0) {
  stop("The model did not converge for lambda = 1.")
}

# Extract the estimated parameters and mean function.
beta_hat <- fit$par
mu_hat1 <- exp(cumsum(beta_hat))

# Plot the observations, true mean and estimated mean.
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

lines(time, mu_hat1,
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

# Define the grid of log-lambda values.
lsp <- seq(-1, 5, by = 0.05)

# Initialise vectors for mBIC and effective degrees of freedom.
mbic <- rep(Inf, length(lsp))
edf <- rep(NA_real_, length(lsp))

# Use the initial parameter estimates for the first fit.
b0 <- beta_start

# Evaluate the model for each penalty value.
for (i in seq_along(lsp)) {
  
  # Convert the log-lambda value to lambda.
  lambda <- exp(lsp[i])
  
  # Fit the model using the current penalty value.
  fit_lambda <- fit_poisson_lasso(y, lambda, b0)
  
  # Skip penalty values where optimisation fails.
  if (fit_lambda$convergence != 0) {
    next
  }
  
  beta_lambda <- fit_lambda$par
  
  edf[i] <- poisson_lasso_edf(beta_lambda, y, lambda)
  
  # Calculate the fitted Poisson mean function.
  mu_lambda <- exp(cumsum(beta_lambda))
  
  # Calculate the maximised log-likelihood.
  ll <- sum(dpois(y, mu_lambda, log = TRUE))
  
  # Calculate the modified Bayesian Information Criterion.
  mbic[i] <- 3 * (edf[i] - 1) * log(n) - 2 * ll
  
  # Use the current estimates as starting values for the next fit.
  b0 <- beta_lambda
}

# Check that at least one penalty value has a valid mBIC.
if (all(!is.finite(mbic))) {
  stop("No valid mBIC values were obtained.")
}

# Find the penalty value that minimises mBIC.
i_opt <- which.min(mbic)
lambda_best <- exp(lsp[i_opt])

# Display the optimal penalty value.
lambda_best

# Display the optimal log-lambda, EDF and mBIC.
data.frame(
  log_lambda = lsp[i_opt],
  lambda = lambda_best,
  EDF = edf[i_opt],
  mBIC = mbic[i_opt]
)
# mBIC is minimised at log(lambda) = 3, giving lambda = 20.08554.

# Plot mBIC against log-lambda.
plot(lsp, mbic,
     type = "l",
     xlab = "log(lambda)",
     ylab = "mBIC",
     main = "Choosing lambda using mBIC")

points(lsp[i_opt], mbic[i_opt],
       col = "red",
       pch = 19)


# Task 5: Fit the final model

# Refit the model using the optimal penalty value.
fit_best <- fit_poisson_lasso(
  y, lambda_best,
  c(log(mean(y)), rep(0, n - 1))
)

# Check convergence of the final model.
fit_best$convergence

if (fit_best$convergence != 0) {
  stop("The final model did not converge.")
}

# Extract the final parameter estimates and mean function.
beta_best <- fit_best$par
mu_hat_best <- exp(cumsum(beta_best))

# Calculate the effective number of changepoints.
edf_best <- poisson_lasso_edf(beta_best, y, lambda_best)
m_hat <- floor(edf_best - 1)

# Identify the changepoints using the largest absolute coefficients.
# Exclude the first parameter because it represents the intercept.
if (m_hat > 0) {
  tau_hat <- sort(order(abs(beta_best[-1]),
                        decreasing = TRUE)[seq_len(m_hat)] + 1)
} else {
  tau_hat <- numeric(0)
}

# Plot the observations, true mean and estimated mean.
plot(time, y,
     pch = 16,
     col = "grey",
     xlab = "Time",
     ylab = "Count / Mean",
     main = "Final Poisson Changepoint Model")

lines(time, mu_true,
      col = "black",
      lwd = 2,
      type = "s")

lines(time, mu_hat_best,
      col = "red",
      lwd = 2,
      type = "s")

abline(v = tau_hat,
       col = "blue",
       lty = 2)

legend("topright",
       cex = 0.55,
       legend = c("Data", "True mean",
                  "Estimated mean", "Changepoints"),
       col = c("grey", "black", "red", "blue"),
       pch = c(16, NA, NA, NA),
       lty = c(NA, 1, 1, 2),
       bty = "n")

# Display the final model results.
lambda_best
m_hat
tau_hat

# EDF = 5.864, giving 4 estimated changepoints.
# Estimated locations: 30, 69, 70 and 72.
# True locations: 30 and 70.
# The change near 70 is split into several estimated changepoints.
# The penalty shrinks the estimated mean changes.


# Task 6: Parametric bootstrap

# Bootstrap sampling is used to estimate uncertainty and bias
# in the fitted Poisson mean function.

# Define the number of bootstrap samples.
B <- 200

# Store the estimated mean function for each bootstrap sample.
mu_boot <- matrix(NA_real_, nrow = B, ncol = n)

# Store the convergence status of each bootstrap fit.
conv_boot <- numeric(B)

# Generate and fit the bootstrap samples.
set.seed(1)

for (j in 1:B) {
  
  # Simulate new observations using the best fitted mean function.
  y_sim <- rpois(n, lambda = mu_hat_best)
  
  # Initialise the parameters for the bootstrap model.
  b_start <- c(log(mean(y_sim)), rep(0, n - 1))
  
  # Fit the model using the optimal penalty value.
  fit_sim <- fit_poisson_lasso(
    y_sim,
    lambda_best,
    beta_start = b_start
  )
  
  # Record the convergence status.
  conv_boot[j] <- fit_sim$convergence
  
  # Skip bootstrap samples where optimisation fails.
  if (fit_sim$convergence != 0) {
    warning(paste("No convergence in bootstrap sample", j))
    next
  }
  
  # Store the estimated mean function.
  mu_boot[j, ] <- exp(cumsum(fit_sim$par))
}

# Check the number of unsuccessful bootstrap fits.
sum(conv_boot != 0)

# Stop if any bootstrap sample failed.
if (any(conv_boot != 0) || anyNA(mu_boot)) {
  stop("Some bootstrap fits failed. Check convergence before continuing.")
}

# Calculate the average bootstrap mean function.
mu_bar <- colMeans(mu_boot)

# Plot the observations, true mean, estimated mean and bootstrap average.
plot(time, y,
     pch = 16,
     col = "grey",
     ylim = c(0, 1.3 * max(y, mu_true, mu_hat_best, mu_bar)),
     xlab = "Time",
     ylab = "Count / Mean",
     main = "Bootstrap Average vs Estimated Mean")

lines(time, mu_true,
      col = "black",
      lwd = 2,
      type = "s")

lines(time, mu_hat_best,
      col = "red",
      lwd = 2,
      type = "s")

lines(time, mu_bar,
      col = "darkgreen",
      lwd = 2,
      type = "s")

legend("topright",
       legend = c("Observations", "True mean",
                  "Estimated mean", "Bootstrap average"),
       col = c("grey", "black", "red", "darkgreen"),
       pch = c(16, NA, NA, NA),
       lty = c(NA, 1, 1, 1),
       lwd = c(NA, 2, 2, 2))

# Calculate the bias-corrected mean function.
mu_bc <- 2 * mu_hat_best - mu_bar

# Calculate the lower and upper 95% bootstrap quantiles.
q_L <- apply(mu_boot, 2, quantile, probs = 0.025)
q_U <- apply(mu_boot, 2, quantile, probs = 0.975)

# Reflect the bootstrap quantiles to obtain bias-corrected limits.
q_L_bc <- pmax(2 * mu_hat_best - q_U, 0)
q_U_bc <- 2 * mu_hat_best - q_L


# Task 7: Final plot

# Calculate suitable vertical axis limits using all plotted values.
ylim <- range(c(y, mu_true, mu_hat_best, mu_bc,
                q_L_bc, q_U_bc), na.rm = TRUE)

# Ensure that the vertical axis starts at zero.
ylim[1] <- min(0, ylim[1])

# Create the plot without drawing observations initially.
plot(time, mu_hat_best,
     type = "n",
     ylim = ylim,
     xlab = "Time k",
     ylab = expression("Count " * y[k] * " / mean " * mu[k]),
     main = "Bias-corrected Poisson Model")

# Draw the area between the bias-corrected confidence limits.
polygon(c(time, rev(time)),
        c(q_L_bc, rev(q_U_bc)),
        col = adjustcolor("blue", alpha.f = 0.15),
        border = NA)

# Add the original observations.
points(time, y,
       pch = 16,
       col = "grey")

# Add the bias-corrected confidence limits.
lines(time, q_L_bc,
      col = "blue",
      lty = 2,
      lwd = 1.5,
      type = "s")

lines(time, q_U_bc,
      col = "blue",
      lty = 2,
      lwd = 1.5,
      type = "s")

# Add the true, estimated and bias-corrected mean functions.
lines(time, mu_true,
      col = "black",
      lwd = 2,
      type = "s")

lines(time, mu_hat_best,
      col = "red",
      lwd = 2,
      type = "s")

lines(time, mu_bc,
      col = "blue",
      lwd = 2,
      type = "s")

# Identify all elements of the final plot.
legend("topright",
       cex = 0.7,
       legend = c("Observations", "True mean",
                  "Estimated mean", "Bias-corrected mean",
                  "Bias-corrected 95% CI"),
       col = c("grey", "black", "red", "blue", "blue"),
       pch = c(16, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 1, 2),
       lwd = c(NA, 2, 2, 2, 1.5),
       bg = "white")

# Interpretation of the final plot:
# Shrinkage pulls the segment means of mu_hat towards each other; the bias
# correction largely reverses this, especially in the middle segment.
round(c(
  mu_hat = mean(mu_hat_best[31:68]),
  mu_bc = mean(mu_bc[31:68])
), 1) 

# The band is narrow within segments but widens sharply around k = 30 and k = 70:
# the bootstrap fits disagree on where the jump is, so mu_bar is smeared and
# mu_bc = 2*mu_hat - mu_bar overshoots just before the jump.
round(mu_bc[29], 1)

# Calculate the proportion of time points where the confidence interval
# contains the true mean function
mean(q_L_bc <= mu_true & mu_true <= q_U_bc)
# Identify time points where the true mean falls outside the confidence interval
which(!(q_L_bc <= mu_true & mu_true <= q_U_bc))
# 0.97 is close to the nominal 95%. The misses are at the jump near k = 70,
# which the fit spreads over k = 69-72.


