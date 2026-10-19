# =============================================================================
# Coin-Toss Pattern Duel
#
# Authors:
# Raneem Alluhaydan S3044929
# Fagan Ganbarli
# Pierrick de Boer S3034691
#
# Overview:
# This script simulates duels between two binary patterns in a sequence
# of coin tosses. It identifies which pattern appears first and estimates
# the probability that player A wins through repeated simulations.
#
# Main functions:
# 1. firsthit(): Find the first occurrence of a pattern.
# 2. sim_duels(): Simulate multiple duels in batches.
# 3. duel(): Simulate a single duel.
# 4. pwin(): Estimate player A's winning probability.
# =============================================================================


# =============================================================================
# 1. Find the first occurrence of a pattern
# =============================================================================

firsthit <- function(x, a) {
  # Find the ending position of the first occurrence of pattern a
  # in each row of x.
  #
  # Convert a vector into a matrix if necessary. Compare the pattern
  # against every possible starting position and record the positions
  # where all elements match. Return the ending position of the first
  # match, or NA if the pattern does not occur.
  #
  # Args:
  #   x: A vector or matrix containing coin tosses.
  #   a: The binary pattern to search for.
  #
  # Returns:
  #   The ending position of the first occurrence in each row,
  #   or NA if no occurrence is found.
  
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }
  
  L <- ncol(x)
  k <- length(a)
  
  # A pattern longer than the sequence cannot occur.
  if (k > L) {
    return(rep(NA, nrow(x)))
  }
  
  # Check every possible starting position.
  hit <- matrix(TRUE, nrow = nrow(x), ncol = L - k + 1)
  
  for (j in seq_along(a)) {
    hit <- hit & x[, j:(L - k + j), drop = FALSE] == a[j]
  }
  
  # Find the first matching position in each row.
  start <- max.col(hit, ties.method = "first")
  start[rowSums(hit) == 0] <- NA
  
  # Convert starting positions into ending positions.
  return(start + k - 1)
}


# =============================================================================
# 2. Simulate multiple duels
# =============================================================================

sim_duels <- function(a, b, p = 0.5, nreps = 1) {
  # Simulate independent duels between patterns a and b.
  #
  # Generate tosses in blocks and check which pattern appears first.
  # Retain the final k - 1 tosses of unfinished sequences between
  # blocks so that patterns crossing block boundaries are detected.
  # Store the outcomes of completed duels and continue with the
  # remaining simulations.
  #
  # Args:
  #   a: Binary pattern representing player A.
  #   b: Binary pattern representing player B.
  #   p: Probability of obtaining a 1 on each toss.
  #   nreps: Number of independent duels to simulate.
  #
  # Returns:
  #   A list containing the winner (1 for A, 0 for B) and the
  #   total number of tosses for each duel.
  
  k <- length(a)
  
  if (k != length(b)) {
    stop("this is not a fair game!")
  }
  
  if (all(a == b)) {
    stop("a and b are identical")
  }
  
  # Initialize the block length and simulation results.
  time <- 0
  B <- 3 * k
  
  active <- 1:nreps
  winner <- rep(NA, nreps)
  tosses <- rep(NA, nreps)
  
  while (length(active) > 0) {
    
    # Generate the initial block of tosses.
    if (time == 0) {
      W <- matrix(
        rbinom(nreps * B, size = 1, prob = p),
        nrow = nreps
      )
      
    } else {
      # Retain the last k - 1 tosses from each active sequence.
      keep <- W[
        ,
        B - (k - 1) + seq_len(k - 1),
        drop = FALSE
      ]
      
      # Generate new tosses for the unfinished duels.
      W_new <- matrix(
        rbinom(
          nrow(keep) * (B - k + 1),
          size = 1,
          prob = p
        ),
        nrow = nrow(keep)
      )
      
      # Append the new tosses to the retained sequence.
      W <- cbind(keep, W_new)
    }
    
    # Find the first occurrence of each player's pattern.
    time_a <- firsthit(W, a)
    time_b <- firsthit(W, b)
    
    # Identify the winner based on the earliest occurrence.
    a_wins <- !is.na(time_a) &
      (is.na(time_b) | time_a < time_b)
    
    b_wins <- !is.na(time_b) &
      (is.na(time_a) | time_b < time_a)
    
    # Store the results of completed duels.
    decided <- a_wins | b_wins
    
    hit_time <- pmin(time_a, time_b, na.rm = TRUE)
    
    winner[active[a_wins]] <- 1
    winner[active[b_wins]] <- 0
    
    tosses[active[decided]] <- time + hit_time[decided]
    
    # Remove completed duels from the active simulations.
    active <- active[!decided]
    W <- W[!decided, , drop = FALSE]
    
    # Advance the number of newly generated tosses.
    time <- time + B - (k - 1)
  }
  
  return(list(winner = winner, tosses = tosses))
}


duel <- function(a,b,p=0.5){
  
  list = sim_duels(a,b,p,nreps=1)
  
  return(c(list[[1]], list[[2]]))
}



pwin <- function(a,b,p = 0.5,nreps = 10000){
  
  list <- sim_duels(a,b,p,nreps)
  
  p_a_wins <- sum(list[[1]]) / length(list[[1]])
  
  return(p_a_wins)
}
  
  


