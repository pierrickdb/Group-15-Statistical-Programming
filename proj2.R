# =============================================================================
# Coin-Toss Pattern Duel
#
# Authors:
# Pierrick de Boer: S3034691
# Fagan Ganbarli: S3015664
# Raneem Alluhaydan: S3044929


# Contributions:
# Pierrick de Boer ( 33%): Task 1-3
# Fagan Ganbarli ( 33%): Task 4-5
# Raneem Alluhaydan ( 33%): Task 6-7
#
# Overview:
# This script simulates duels between two binary patterns in a sequence
# of coin tosses. It identifies which pattern appears first and estimates
# the probability that player A wins through repeated simulations.
#
# Main functions:
# 1. firsthit(): Find the first occurrence of a pattern.
# 2. sim_duels(): Simulate multiple duels.
# 3. duel(): Simulate one duel.
# 4. pwin(): Estimate the probability that a pattern wins.
# 5. ewait(): Calculate the expected waiting time of a pattern.
# 6. pmat(): Calculate win probabilities between patterns.
# 7. find_cycle(): Find a non-transitive cycle of patterns.
# =============================================================================

# Task 1: Find the first occurrence of a pattern in a toss stream.

# Idea:
# Compare the pattern with every window of k consecutive tosses at
# once and return the position where the first matching window ends.

# Input:
# x = A toss stream, either a vector or a matrix with one
#     stream per row.
# a = A pattern of length k.

# Output:
# The toss at which the pattern is first completed in each row,
# or NA if the pattern does not appear.

# Method:
# Convert a vector into a one-row matrix.
# Build a logical matrix of all possible starting positions and
# compare each element of the pattern with the matching shifted
# columns of x.
# Take the first matching start in each row and add k - 1 to
# obtain the ending position.

firsthit <- function(x, a) {
  
  # Treat a single stream as a matrix with one row.
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

# Test firsthit() function.
x <- c(1, 0, 0, 1, 0, 1, 1, 1, 0)

firsthit(x, c(1, 0, 1))   # Expected: 6
firsthit(x, c(0, 1, 1))   # Expected: 7
firsthit(x, c(0, 0, 0))   # Expected: NA


# Task 2: Simulate duels between two patterns.

# Idea:
# Toss a shared coin stream for many duels at once and stop each
# duel as soon as one of the two patterns appears.

# Input:
# a = Binary pattern of player A.
# b = Binary pattern of player B, distinct from a and of equal length.
# p = The probability of a head
# nreps = The number of independent duels to simulate.

# Output:
# A list containing the winner of each duel (1 for A, 0 for B)
# and the number of tosses each duel took.

# Method:
# Check that the patterns have equal length and are distinct.
# Generate tosses in blocks, one row per duel, and find the first
# occurrence of each pattern with firsthit().
# Store the results of decided duels and continue only with the
# undecided ones. Keep the last k - 1 tosses of each undecided row
# so that patterns crossing block boundaries are detected.

sim_duels <- function(a, b, p = 0.5, nreps = 1) {
  
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
      keep <- W[,B - (k - 1) + seq_len(k - 1),drop = FALSE]
      
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


# Task 2: Simulate a single duel.

# Idea:
# A single duel is a special case of sim_duels() with one repetition,
# so the simulation code is shared rather than repeated.

# Input:
# a = Binary pattern of player A.
# b = Binary pattern of player B.
# p = The probability of a head

# Output:
# A vector c(winner, tosses), where winner is 1 if A wins and 0
# if B wins, and tosses is the number of tosses the duel took.

# Method:
# Call sim_duels() with nreps = 1 and combine both results
# into one vector.

duel <- function(a, b, p = 0.5) {
  
  result <- sim_duels(a, b, p, nreps = 1)
  
  return(c(result[[1]], result[[2]]))
}


# Task 3: Estimate the probability that one pattern wins a duel.

# Idea:
# Estimate the win probability of player A as the proportion of
# simulated duels that A wins.

# Input:
# a = Binary pattern of player A.
# b = Binary pattern of player B.
# p = The probability of a head (0 < p < 1).
# nreps = The number of independent duels to simulate.

# Output:
# The estimated probability that pattern a appears before pattern b.

# Method:
# Simulate all duels at once with sim_duels().
# Since a win for A is stored as 1 and a loss as 0, the mean of
# the winner vector estimates the win probability.

pwin <- function(a, b, p = 0.5, nreps = 10000) {
  
  # Simulate the duels.
  result <- sim_duels(a, b, p, nreps)
  
  # Calculate the proportion of duels won by player A.
  p_a_wins <- sum(result[[1]]) / length(result[[1]])
  
  return(p_a_wins)
}


# Task 4.a: Calculating the exact expected waiting time for a pattern.

# Idea:
# The expected waiting time depends on the lengths for which
# the first and last j elements of a pattern match. Each matching
# overlap contributes the reciprocal of the probability of its prefix.

# Input:
# a = A non-empty binary vector representing the pattern.
# p = The probability of a head (0 < p < 1).

# Output:
# The exact expected waiting time until the pattern appears. 

# Method:
# Check that the inputs are valid.
# Examine every possible overlap length. If the first and last j
# elements match, calculate the prefix probability and add its
# reciprocal to the expected waiting time.

ewait <- function(a, p = 0.5) {
  
  # Check that the pattern is a non-empty binary vector.
  if (length(a) == 0 || any(a != 0 & a != 1) || anyNA(a)) {
    stop("The pattern must be a non-empty binary vector.")  
  }
  
  # Check that p is a valid probability between 0 and 1.
  if (!is.numeric(p) || length(p) != 1 || 
      p <= 0 || p >= 1 || is.na(p)) {
    stop("p must be a numeric value between 0 and 1.")  
  }
  
  k <- length(a)
  expected_wait <- 0
  
  # Examine every possible overlap length, including the full pattern.
  for (j in 1:k) {
    
    # An overlap exists when the first and last j elements match.
    if (all(a[1:j] == a[(k-j+1):k])) {
      
      # Count the number of heads in the prefix.
      n_head <- sum(a[1:j] == 1)
      
      # Calculate the probability of the prefix.
      r_probability <- (p ^ n_head) * ((1 - p) ^ (j - n_head))
      
      # Add the reciprocal probability for this overlap.
      expected_wait <- expected_wait + (1 / r_probability)
    }
  }
  return(expected_wait)
}

# Check the expected waiting times for a fair coin.
patterns <- list(
  HH = c(1, 1),
  TT = c(0, 0),
  HT = c(1, 0),
  TH = c(0, 1)
)

exact <- sapply(patterns, ewait)
print(exact)

# Task 4.b: Calculate the expected waiting time using Monte Carlo simulation.

# Idea:
# The function simulates coin tosses until the given pattern appears.
# It repeats this process nreps times and returns the average.

# Input:
# a = A binary pattern.
# p = Probability of getting a head.
# nreps = Number of simulations.

# Output:
# Monte-Carlo estimated expected waiting time.

# Method:
# Generate tosses and use firsthit() to find the pattern.
# Calculate the average of results.
# Compare them with the exact expected waiting time.

mc_ewait <- function(a, p = 0.5, nreps = 10000) {
  
  waiting_times <- numeric(nreps)
  
  # Run simulation nreps times. 
  for (i in seq_len(nreps)) {
    
    # Create an empty vector.
    stream <- integer(0)
    
    repeat {
      
      # Add one toss to the stream.
      stream <- c(stream, rbinom(1, 1, p))
      
      # Check whether the pattern has appeard.
      hit <- firsthit(stream, a)
      
      # If pattern has appeard, stop the repeat.
      if (!is.na(hit[1])) {
        waiting_times[i] <- hit[1]
        break
      }
    }
  }
  
  return(mean(waiting_times))
}

# Compare Monte Carlo estimates with the exact waiting values 
# for patterns of length 2 with a fair coin.
set.seed(7)
mc <- sapply(patterns,  mc_ewait, p = 0.5, nreps = 10000)

results <- data.frame(
  Pattern = names(patterns),
  Exact = as.numeric(exact),
  Monte_Carlo = as.numeric(mc)
)

print(results)

# The results show that the Monte Carlo estimates are close to the 
# exact expected waiting times for all four patterns. For a fair 
# coin, the expected waiting times are 6 tosses for HH and TT, 
# and 4 tosses for HT and TH. 

# Compare Monte Carlo estimates with the exact waiting values 
# for a biased coin: p is not equal to 0.5.
exact_biased <- ewait(c(1, 1, 0), p = 0.7)
mc_biased <- mc_ewait(c(1, 1, 0), p = 0.7, nreps = 10000)

Biased_coin <- data.frame(
  Biased_Probability = 0.70,
  Exact_Waiting = exact_biased,
  Monte_Carlo_Waiting = mc_biased
)

print(Biased_coin)

# For a biased coin with p = 0.7, the exact expected waiting time 
# for HHT is approximately 6.8027, while the Monte Carlo estimate 
# is 6.8034. These results confirm that the simulation provides 
# a good approximation of the exact expected waiting time.


# Task 5: Calculate the matrix of pairwise win probabilities.

# Idea:
# Estimate how often each pattern beats every other pattern using 
# repeated simulations, and store the result probabilities in a matrix. 

# Input:
# k = The length of each binary pattern.
# p = The probability of a head (0 < p < 1).
# nreps = The number of simulated duels for each pair..

# Output:
# A 2^k by 2^k matrix of estimated win probabilities.
# The diagonal entries are NA, because a pattern cannot duel itself.

# Method:
# Generate all possible patterns of length k and arrange them in
# lexicographic order, with 0 before 1.
# Estimate the win probability for every distinct pattern using 
# the pwin() function and store the result in matrix M.

pmat <- function(k, p = 0.5, nreps = 10000) {
  
  # Generate all possible binary patterns of length k.
  pattern <- expand.grid(rep(list(0:1), k))
  
  # Arrange patterns in lexicographic order, with 0 before 1.
  pattern <- pattern[
    do.call(order, as.data.frame(pattern)), , drop = FALSE]
  
  # Convert the patterns to a matrix.
  pattern <- as.matrix(pattern)
  
  # Create a matrix to store the estimated win probabilities.
  M <- matrix(NA, nrow = 2^k, ncol = 2^k)
  
  # Calculate the number of patterns.
  n <- nrow(pattern)
  
  # Label the rows and columns using the binary patterns.
  rownames(M) <- apply(pattern, 1, paste, collapse = "")
  colnames(M) <- apply(pattern, 1, paste, collapse = "")
  
  # Estimate the probability that each pattern beats every other pattern.
  for (i in 1:n) {
    for (j in 1:n) {
      
      # Leave diagonal entries as NA, because a pattern cannot duel itself.
      if (i != j) {
        M[i, j] <- pwin(
          pattern[i, ],
          pattern[j, ],
          p = p,
          nreps = nreps
        )
      }
    }
  }
  return(M)
}



# Task 5(a): Find the best reply and its win probability for each pattern.

# Idea:
# For each pattern, find the pattern that beats it with the 
# highest probability and display the results in a table.

# Input:
# M = The matrix of win probabilities created by pmat() function
# with k = 3 and a fair coin.

# Output:
# A data frame showing each pattern, its best reply, and the 
# probability that the best reply wins.

# Method:
# Replace the diagonal NA values with -inf so they cannot be 
# selected as the best reply.
# Find the maximum value in each column.
# Record the each pattern, its best reply, and the associated 
# win probability in a data frame.

# Create a matrix for k = 3 and a fair coin using pmat() function.
M <- pmat(k = 3, p = 0.50)

# Replace the diagonal NA values with -inf.
M <- replace(M, is.na(M), -Inf)

# Find the index of the best reply for each pattern.
best_reply <- apply(M, 2, which.max)

# Create a vector to store the win probabilities.
win_probability <- numeric(ncol(M))

# Store the win probability of each best reply.
for (j in seq_len(ncol(M))) {
  win_probability[j] <- M[best_reply[j], j]
}

# Create a data frame of patterns, best replies, and win probabilities.
best_reply_table <- data.frame(
  Pattern = rownames(M),
  Best_reply = rownames(M)[best_reply],
  Probability = win_probability
)

print(best_reply_table)


# Task 5(b)(1): Find evenly matched pattern pairs.

# Idea:
# Find pairs with approximately equal chances (between 0.49 and 0.51) 
# of winning against each other.

# Input:
# M = The matrix of win probabilities created by pmat() function
# with k = 3 and a fair coin.

# Output:
# A data frame showing each evenly matched pair and the win 
# probabilities of both pattern

# Method:
# Examine each distinct pair only once.
# Keep pairs for which both win probabilities are between 0.49 and 0.51.

# Create an empty data frame to store the results.
evenly_matched_pairs <- data.frame(
  Pattern1 = character(),
  Pattern2 = character(),
  P1_Win_Probability = numeric(),
  P2_Win_Probability = numeric()
)

# Examine all pairs with for loop.
for (i in seq_len(nrow(M))) {
  for (j in seq_len(ncol(M))) {
    
    # Check each pair once and exclude the diagonal.
    if (j > i && 
        !is.na(M[i, j]) &&
        !is.na(M[j, i]) &&
        abs(M[i, j] - 0.5) <= 0.01 &&
        abs(M[j, i] - 0.5) <= 0.01) {
      
      # Store both patterns and their win probabilities.
      evenly_matched_pairs <- rbind(
        evenly_matched_pairs,
        data.frame(
          Pattern1 = rownames(M)[i],
          Pattern2 = colnames(M)[j],
          P1_Win_Probability = M[i, j],
          P2_Win_Probability = M[j, i]
        )
      )
    }
  }
}

# Reset row names and display the result
row.names(evenly_matched_pairs) <- NULL
print(evenly_matched_pairs)


# Task 5(b)(2): Find the shortest non-transitive cycle.

# Idea:
# Find the shortest cycle in which every pattern beats the next
# pattern with at least given win probability.

# Input:
# m = A matrix of pairwise win probabilities.
# min_probability = The minimum probability required for each win.

# Output:
# A data frame showing the patterns in each cycle, the direction 
# of each win, and win probabilities.

# Method:
# Search cycles in increasing order of length. Use recursive search
# to extend path without repeating patterns, and check whether the
# final pattern beats the starting pattern to complete a cycle.

find_cycle <- function(m, min_probability = 0.60) {
  
  k <- nrow(m)
  
  # Check cycle from the shortest possible length (3).
  for (i in 3:k) {
    
    cycles <- list()
    
    # Try each pattern as the starting point.
    for(start in seq_len(k)) {
      
      search <- function(path) {
        
        last <- path[length(path)]
        
        # Check whether the path can be closed into a cycle.
        if (length(path) == i) {
          
          if (!is.na(m[last, start]) &&
              m[last, start] >= min_probability) {
            
            cycles[[length(cycles) + 1]] <<- path
          }
          
          return(invisible(NULL))
        }
        
        # Only consider unused patterns to avoid duplicate cycles.
        sequence <- setdiff(seq_len(k), path)
        sequence <- sequence[sequence > start]
        
        for (next_item in sequence) {
          
          if (!is.na(m[last, next_item]) && 
              m[last, next_item] >= min_probability) {
            
            search(c(path, next_item))
          }
        }
      }
      search(start)
    }
    
    # Return the cycle when shortest length has been found.
    if (length(cycles) > 0) {
      
      result <- do.call(rbind, lapply(seq_along(cycles), function(x) {
        
        path <- cycles[[x]]
        next_path <- c(path[-1], path[1])
        
        data.frame(
          Cycle = x,
          From = rownames(m)[path],
          To = colnames(m)[next_path],
          Win = m[cbind(path, next_path)]
        )
      }))
      rownames(result) <- NULL
      return(result)
    }
  }
  
  # message dos not exist 
  return("No cycle found.")
}

# Find the shortest cycle with minimum win probability = 0.60.
find_cycle(M, min_probability = 0.60)  

# Task 5(c): Compare expected waiting times with duel results.

# Idea:
# Find pairs with the same expected waiting time and identify.
# pairs with noticeably different duel probabilities.

# Input:
# M = The matrix of win probabilities created by pmat() function
# with k = 3 and a fair coin.
# ewait() = The function that calculates the expected waiting time. 

# Output:
# A data frame that demonstrates pairs with the same expected 
# waiting time whose duel probability is below 0.4 or above 0.6.

# Method:
# Get the pattern labels from the row names of M.
# Convert labels into a binary matrix.
# Calculate the exact expected waiting time with ewait() function 
# Compare each pairs and keep whose expected waiting times are the same.
# Display pairs with win probability notably far from 0.5.

# Get the pattern labels from the row names of M.
pattern_label <- rownames(M)

# Convert labels into a binary matrix.
pattern_label <- t(sapply(pattern_label, function(x) {
  as.numeric(strsplit(x, "")[[1]])
})) 

# Define pattern names.
pattern_names <- apply(pattern_label, 1, paste, collapse = "")

# Calculate th exact expected waiting time.
expected_wait <- apply(pattern_label, 1, ewait)

# Create an empty data frame to store pairs.
comparison_table <- data.frame(
  Pattern1 = character(),
  Pattern2 = character(),
  Expected_Wait = numeric(),
  Win_Probability = numeric()
)

# Compare every distinct pairs.
for (i in seq_len(nrow(M))) {
  for (j in seq_len(ncol(M))) {
    
    # Check each pair ony one time
    # Keep pairs with the same expected waiting time.
    if (i < j &&
        expected_wait[i] == expected_wait[j]) {
      
      # Store the pair, their equal expected waiting time,
      # and the probability that Pattern1 beats Pattern2.
      comparison_table <- rbind(
        comparison_table, data.frame(
          Pattern1 = pattern_names[i],
          Pattern2 = pattern_names[j],
          Expected_Wait = expected_wait[i],
          Win_Probability = M[i, j]
        )
      )
    }
  }
}

# Display pairs whose win probabilities notable differ from 0.5 
# (win probability is below 0.4 or above 0.6.)

comparison_table <- comparison_table[
  abs(comparison_table$Win_Probability - 0.5) > 0.1, ]

row.names(comparison_table) <- NULL

print(comparison_table)

# Task 6: Plot duel lengths and variation in win probabilities.

# Idea: Investigate duel lengths and how win-probability estimates vary.
# Input: Patterns 011 and 111, a fair coin, and 2000 duels per run.
# Output: A graph of duel lengths and a graph of 50 win-probability estimates.
# Method: Simulate duels, mark the median length, repeat the win-probability
# estimate 50 times, and mark the mean estimate.

# Simulate 2000 duels between patterns 011 and 111.
a <- c(0, 1, 1)
b <- c(1, 1, 1)
results <- sim_duels(a, b, p = 0.5, nreps = 2000)

# Plot the distribution of duel lengths and mark the median.
tosses <- results$tosses
plot(ecdf(tosses), verticals = TRUE, do.points = FALSE,
     main = "Distribution of duel lengths",
     xlab = "Number of tosses",
     ylab = "Estimated probability")

med <- median(tosses)
lines(c(med, med), c(0, 1))
text(med + 2, 0.95, "Median")

# Estimate the win probability in 50 independent runs.
win_estimates <- numeric(50)

for (i in 1:50) {
  win_estimates[i] <- pwin(
    a = c(0, 1, 1),
    b = c(1, 1, 1),
    p = 0.5,
    nreps = 2000
  )
}

# Calculate the mean estimate and plot the results.
mean_estimate <- mean(win_estimates)

plot(win_estimates,
     xlab = "Run",
     ylab = "Win probability",
     main = "Win probability across 50 runs")

abline(h = mean_estimate, lty = 2)
text(43, mean_estimate, "Mean", pos = 3)


# Task 7: Investigate the effect of coin bias.

# Idea: Check how changing the coin bias affects the best replies.
# Input: Patterns of length 3, with p = 0.25 and p = 0.7.
# Output: A best-reply table for each coin bias.
# Method: Calculate the win-probability matrix, find the best reply
# to each pattern, and compare the results for both biases.

# Analyse the best replies for p = 0.25.
M <- pmat(k = 3, p = 0.25, nreps = 2000)
M <- replace(M, is.na(M), -Inf)

best_reply <- apply(M, 2, which.max)
win_probability <- numeric(ncol(M))

for (j in seq_len(ncol(M))) {
  win_probability[j] <- M[best_reply[j], j]
}

best_reply_table <- data.frame(
  Pattern = colnames(M),
  Best_reply = rownames(M)[best_reply],
  Probability = win_probability,
  row.names = NULL
)

print(best_reply_table)


# Repeat the analysis for p = 0.7.
M <- pmat(k = 3, p = 0.7, nreps = 2000)
M <- replace(M, is.na(M), -Inf)

best_reply <- apply(M, 2, which.max)
win_probability <- numeric(ncol(M))

for (j in seq_len(ncol(M))) {
  win_probability[j] <- M[best_reply[j], j]
}

best_reply_table <- data.frame(
  Pattern = colnames(M),
  Best_reply = rownames(M)[best_reply],
  Probability = win_probability,
  row.names = NULL
)

print(best_reply_table)


# Best replies and cycles:
# For p = 0.25, one observed cycle is:
# 100 -> 000 -> 001 -> 010 -> 100
#
# For p = 0.7, one observed cycle is:
# 101 -> 011 -> 111 -> 110 -> 101
#
# Changing the coin bias changes the best replies and the patterns
# involved in the observed cycles.

# Overall analysis:
# Ranking patterns by their individual waiting times does not
# necessarily predict who wins a duel.
# Both patterns are observed in the same sequence of coin tosses,
# so their appearances are not independent.
# Therefore, isolated waiting-time summaries can be misleading
# when predicting duel outcomes.
