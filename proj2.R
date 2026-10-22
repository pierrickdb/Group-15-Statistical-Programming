# =============================================================================
# Coin-Toss Pattern Duel
#
# Authors:
# Pierrick de Boer: S3034691
# Fagan Ganbarli: S3015664
# Raneem Alluhaydan: S3044929
#
# Contributions:
# Pierrick de Boer (33%): Task 1-3
# Fagan Ganbarli (33%): Task 4-5
# Raneem Alluhaydan (33%): Task 6-7
#
# Overview:
# This script simulates duels between two binary patterns in a sequence
# of coin tosses. It identifies which pattern appears first and estimates
# the probability that player A wins through repeated simulations.
# Heads are coded as 1 and tails as 0.
#
# Main functions:
# 1. firsthit(): Find the first occurrence of a pattern.
# 2. sim_duels(): Simulate multiple duels.
# 3. duel(): Simulate one duel.
# 4. pwin(): Estimate the probability that a pattern wins.
# 5. ewait(): Calculate the exact expected waiting time of a pattern.
# 6. mc_ewait(): Estimate the expected waiting time by simulation.
# 7. pmat(): Calculate win probabilities between all patterns.
# 8. find_cycle(): Find a non-transitive cycle of patterns.
# =============================================================================


# Task 1: Find the first occurrence of a pattern in a toss stream.

firsthit <- function(x, a) {
  # Finds the toss at which pattern a is completed for the first time
  # in each toss stream.
  #
  # Inputs:
  #   x: toss stream, either a vector or a matrix with one stream per row
  #   a: binary pattern of length k
  #
  # Output:
  #   The ending position of the first occurrence in each row,
  #   or NA if the pattern does not appear.
  #
  # Method:
  #   Convert a vector into a one-row matrix. Build a logical matrix
  #   of all possible starting positions and compare each element of
  #   the pattern with the matching shifted columns of x. Take the
  #   first matching start in each row and add k - 1 to obtain the
  #   ending position.
  
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

# Test firsthit() on the stream H T T H T H H H T, which can be
# checked by hand: HTH ends at toss 6, THH at toss 7, TTT never occurs.
x <- c(1, 0, 0, 1, 0, 1, 1, 1, 0)

firsthit(x, c(1, 0, 1))   # Expected: 6
firsthit(x, c(0, 1, 1))   # Expected: 7
firsthit(x, c(0, 0, 0))   # Expected: NA


# Task 2: Simulate duels between two patterns.

sim_duels <- function(a, b, p = 0.5, nreps = 1) {
  # Simulates independent duels in which both players watch the same
  # coin stream, and stops each duel as soon as one pattern appears.
  #
  # Inputs:
  #   a: binary pattern of player A
  #   b: binary pattern of player B, distinct from a and of equal length
  #   p: probability of a head (0 < p < 1)
  #   nreps: number of independent duels to simulate
  #
  # Output:
  #   A list containing the winner of each duel (1 for A, 0 for B)
  #   and the number of tosses each duel took.
  #
  # Method:
  #   Check that the patterns have equal length and are distinct.
  #   Generate tosses in blocks, one row per duel, and find the first
  #   occurrence of each pattern with firsthit(). Store the results of
  #   decided duels and continue only with the undecided ones. Keep the
  #   last k - 1 tosses of each undecided row so that patterns crossing
  #   block boundaries are detected.
  
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

duel <- function(a, b, p = 0.5) {
  # Simulates one duel. A single duel is a special case of sim_duels()
  # with one repetition, so the simulation code is shared.
  #
  # Inputs:
  #   a: binary pattern of player A
  #   b: binary pattern of player B
  #   p: probability of a head (0 < p < 1)
  #
  # Output:
  #   A vector c(winner, tosses), where winner is 1 if A wins and 0
  #   if B wins, and tosses is the number of tosses the duel took.
  #
  # Method:
  #   Call sim_duels() with nreps = 1 and combine both results
  #   into one vector.
  
  result <- sim_duels(a, b, p, nreps = 1)
  
  return(c(result[[1]], result[[2]]))
}


# Task 3: Estimate the probability that one pattern wins a duel.

pwin <- function(a, b, p = 0.5, nreps = 10000) {
  # Estimates the probability that pattern a appears before pattern b.
  #
  # Inputs:
  #   a: binary pattern of player A
  #   b: binary pattern of player B
  #   p: probability of a head (0 < p < 1)
  #   nreps: number of independent duels to simulate
  #
  # Output:
  #   The estimated probability that player A wins.
  #
  # Method:
  #   Simulate all duels at once with sim_duels(). Since a win for A
  #   is stored as 1 and a loss as 0, the mean of the winner vector
  #   estimates the win probability.
  
  # Simulate the duels.
  result <- sim_duels(a, b, p, nreps)
  
  # Calculate the proportion of duels won by player A.
  p_a_wins <- sum(result[[1]]) / length(result[[1]])
  
  return(p_a_wins)
}


# Task 4(a): Calculate the exact expected waiting time for a pattern.

ewait <- function(a, p = 0.5) {
  # Calculates the exact expected waiting time of a single pattern
  # from its self-overlap structure.
  #
  # Inputs:
  #   a: non-empty binary vector representing the pattern
  #   p: probability of a head (0 < p < 1)
  #
  # Output:
  #   The exact expected number of tosses until the pattern appears.
  #
  # Method:
  #   Check that the inputs are valid. Examine every overlap length j
  #   from 1 to k. If the first and last j elements match, calculate
  #   the probability of the first j outcomes and add its reciprocal
  #   to the expected waiting time.
  
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

# Check the exact expected waiting times of the length-2 patterns
# for a fair coin (expected: 6 for HH and TT, 4 for HT and TH).
patterns <- list(
  HH = c(1, 1),
  TT = c(0, 0),
  HT = c(1, 0),
  TH = c(0, 1)
)

exact <- sapply(patterns, ewait)
print(exact)


# Task 4(b): Estimate the expected waiting time by simulation.

mc_ewait <- function(a, p = 0.5, nreps = 10000) {
  # Estimates the expected waiting time of a pattern by Monte Carlo
  # simulation, as a check on ewait().
  #
  # Inputs:
  #   a: binary pattern
  #   p: probability of a head (0 < p < 1)
  #   nreps: number of simulated waiting times
  #
  # Output:
  #   The Monte Carlo estimate of the expected waiting time.
  #
  # Method:
  #   For each repetition, add one toss at a time to the stream and
  #   use firsthit() to check whether the pattern has appeared. Record
  #   the waiting time once it appears and return the average over
  #   all repetitions.
  
  waiting_times <- numeric(nreps)
  
  # Run the simulation nreps times.
  for (i in seq_len(nreps)) {
    
    # Start with an empty stream.
    stream <- integer(0)
    
    repeat {
      
      # Add one toss to the stream.
      stream <- c(stream, rbinom(1, 1, p))
      
      # Check whether the pattern has appeared.
      hit <- firsthit(stream, a)
      
      # Stop once the pattern has appeared.
      if (!is.na(hit[1])) {
        waiting_times[i] <- hit[1]
        break
      }
    }
  }
  
  return(mean(waiting_times))
}

# Compare the Monte Carlo estimates with the exact waiting times
# for the length-2 patterns and a fair coin.
set.seed(7)
mc <- sapply(patterns,  mc_ewait, p = 0.5, nreps = 10000)

results <- data.frame(
  Pattern = names(patterns),
  Exact = as.numeric(exact),
  Monte_Carlo = as.numeric(mc)
)

print(results)

# The Monte Carlo estimates are close to the exact expected waiting
# times for all four patterns: 6 tosses for HH and TT, and 4 tosses
# for HT and TH.

# Compare the Monte Carlo estimate with the exact waiting time
# for a biased coin (p = 0.7).
exact_biased <- ewait(c(1, 1, 0), p = 0.7)
mc_biased <- mc_ewait(c(1, 1, 0), p = 0.7, nreps = 10000)

Biased_coin <- data.frame(
  Biased_Probability = 0.70,
  Exact_Waiting = exact_biased,
  Monte_Carlo_Waiting = mc_biased
)

print(Biased_coin)

# For p = 0.7, the exact expected waiting time for HHT is about
# 6.8027 and the Monte Carlo estimate is 6.8034, so the simulation
# confirms the exact result for a biased coin as well.


# Task 5: Calculate the matrix of pairwise win probabilities.

pmat <- function(k, p = 0.5, nreps = 10000) {
  # Estimates how often each pattern of length k beats every other
  # pattern of the same length in a duel.
  #
  # Inputs:
  #   k: length of each binary pattern
  #   p: probability of a head (0 < p < 1)
  #   nreps: number of simulated duels for each pair
  #
  # Output:
  #   A 2^k by 2^k matrix in which entry (i, j) is the estimated
  #   probability that pattern i beats pattern j. The diagonal is NA,
  #   because a pattern cannot duel itself. Rows and columns are
  #   labelled by the patterns, e.g. "010".
  #
  # Method:
  #   Generate all 2^k patterns and sort them in lexicographic order,
  #   with 0 before 1 (000, 001, 010, ..., 111). Estimate the win
  #   probability of every ordered pair of distinct patterns with
  #   pwin() and store it in the matrix.
  
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
#
# For each pattern (column of M), the best reply is the pattern that
# beats it with the highest probability. The diagonal NA values are
# replaced by -Inf so that a pattern cannot be its own best reply.

# Create the win-probability matrix for k = 3 and a fair coin.
M <- pmat(k = 3, p = 0.50)

# Replace the diagonal NA values with -Inf.
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
#
# A pair is evenly matched if both patterns win with a probability
# between 0.49 and 0.51. Each distinct pair is examined only once.

# Create an empty data frame to store the results.
evenly_matched_pairs <- data.frame(
  Pattern1 = character(),
  Pattern2 = character(),
  P1_Win_Probability = numeric(),
  P2_Win_Probability = numeric()
)

# Examine all pairs.
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

# Reset row names and display the result.
row.names(evenly_matched_pairs) <- NULL
print(evenly_matched_pairs)


# Task 5(b)(2): Find the shortest non-transitive cycle.

find_cycle <- function(m, min_probability = 0.60) {
  # Finds the shortest cycle P1 -> P2 -> ... -> Pm -> P1 of distinct
  # patterns in which every pattern beats the next one.
  #
  # Inputs:
  #   m: matrix of pairwise win probabilities, as returned by pmat()
  #   min_probability: minimum win probability required for each step
  #
  # Output:
  #   A data frame listing every shortest cycle found, with the
  #   direction of each win and its probability, or a message if no
  #   cycle exists.
  #
  # Method:
  #   Search cycle lengths in increasing order, starting at 3. For each
  #   starting pattern, extend the path recursively with unused patterns
  #   that the last pattern beats, and check whether the final pattern
  #   beats the starting pattern. Only patterns with a larger index than
  #   the start are used, so each cycle is found only once.
  
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
  
  # Return a message if no cycle exists.
  return("No cycle found.")
}

# Find the shortest cycle with minimum win probability = 0.60.
find_cycle(M, min_probability = 0.60)  


# Task 5(c): Compare expected waiting times with duel results.
#
# Calculate the exact expected waiting time of every pattern with
# ewait(), keep the pairs with equal waiting times, and show those
# whose duel win probability is below 0.4 or above 0.6.

# Get the pattern labels from the row names of M.
pattern_label <- rownames(M)

# Convert labels into a binary matrix.
pattern_label <- t(sapply(pattern_label, function(x) {
  as.numeric(strsplit(x, "")[[1]])
})) 

# Define pattern names.
pattern_names <- apply(pattern_label, 1, paste, collapse = "")

# Calculate the exact expected waiting time.
expected_wait <- apply(pattern_label, 1, ewait)

# Create an empty data frame to store pairs.
comparison_table <- data.frame(
  Pattern1 = character(),
  Pattern2 = character(),
  Expected_Wait = numeric(),
  Win_Probability = numeric()
)

# Compare every distinct pair.
for (i in seq_len(nrow(M))) {
  for (j in seq_len(ncol(M))) {
    
    # Check each pair only once and keep pairs with
    # the same expected waiting time.
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

# Keep pairs whose win probability clearly differs from 0.5
# (below 0.4 or above 0.6).
comparison_table <- comparison_table[
  abs(comparison_table$Win_Probability - 0.5) > 0.1, ]

row.names(comparison_table) <- NULL

print(comparison_table)


# Task 6: Plot duel lengths and variation in win probabilities.
#
# First, simulate 2000 duels between 011 and 111 with a fair coin and
# plot the distribution of duel lengths, marking the median. Second,
# estimate the win probability of 011 (the best reply to 111) in 50
# independent runs of 2000 duels and plot the estimates around their
# mean.

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
#
# Repeat the best-reply analysis of Task 5(a) for p = 0.25 and
# p = 0.7 and compare the results with the fair coin.

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
