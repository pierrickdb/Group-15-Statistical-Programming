# =============================================================================
# Coin-Toss Pattern Duel
#
# Authors:
# Pierrick de Boer S3034691
# Fagan Ganbarli s3015664
# Pierrick de Boer S3034691
# Raneem Alluhaydan S3044929
#
# Contributions:
# Pierrick de Boer ( %): Task 1-3
# Fagan Ganbarli ( %): Task 4-5
# Raneem Alluhaydan ( %): Task 6-7
#
# Overview:
# This script simulates duels between two binary patterns in a sequence
# of coin tosses. It identifies which pattern appears first and estimates
# the probability that player A wins through repeated simulations.
#
# Main functions:
# 1. firsthit(): Find the first occurrence of a pattern. (Pierrick)
# 2. sim_duels(): Simulate multiple duels in batches.(Pierrick)
# 3. duel(): Simulate a single duel. (Pierrick)
# 4. pwin(): Estimate player A's winning probability.(Pierrick)
# 5. 
# 6. 
# 7.
# 8.
# =============================================================================

# Task 1: 

# Idea:
#

# Input:
#

# Output:
#

# Method:
#

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


# Task 2: 

# Idea:
#

# Input:
#

# Output:
#

# Method:
#

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



duel <- function(a, b, p = 0.5) {
  # Simulate one duel by calling sim_duels() with a single repetition.
  #
  # Args:
  #   a: Binary pattern representing player A.
  #   b: Binary pattern representing player B.
  #   p: Probability of obtaining a 1 on each toss.
  #
  # Returns:
  #   A vector containing the winner and the total number of tosses.
  
  result <- sim_duels(a, b, p, nreps = 1)
  
  return(c(result[[1]], result[[2]]))
}


# Task 3: 

# Idea:
#

# Input:
#

# Output:
#

# Method:
#

pwin <- function(a, b, p = 0.5, nreps = 10000) {
  # Estimate player A's probability of winning through simulation.
  #
  # Simulate independent duels and calculate the proportion won
  # by player A. Since wins are represented by 1 and losses by 0,
  # the mean of the winner vector estimates the winning probability.
  #
  # Args:
  #   a: Binary pattern representing player A.
  #   b: Binary pattern representing player B.
  #   p: Probability of obtaining a 1 on each toss.
  #   nreps: Number of independent duels to simulate.
  #
  # Returns:
  #   The estimated probability that player A wins.
  
  result <- sim_duels(a, b, p, nreps)
  
  p_a_wins <- sum(result[[1]]) / length(result[[1]])
  
  return(p_a_wins)
}
  

# Task 4: Calculating the exact expected waiting time for a pattern.

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
  Probablity = win_probability
)

print(best_reply_table)


# Task 5(b)(I): Find evenly matched pattern pairs.

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


# Task 5(b)(II): Find the shortest nontransitive cycle.

# Idea:
# 

# Input:
# 

# Output:
# 

# Method:
# 

find_cycle <- function(m, min_probability = 0.60) {
  
  k <- nrow(m)
  
  for (i in 3:k) {
    
    cycles <- list()
    
    for(start in seq_len(k)) {
      search <- function(path) {
        last <- path[length(path)]
        
        if (length(path) == i) {
          if (!is.na(m[last, start]) &&
              m[last, start] >= min_probability) {
            cycles[[length(cycles) + 1]] <<- path
          }
          return(invisible(NULL))
        }
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
  
  return("No Respound")
}
  

# Task 5(c): Compare expected waiting times with duel results.

# Idea:
# Find pairs with the same expected waiting time and identify.
# pairs with noticeably different duel probabilities.

# Input:
# M = The matrix of win probabilities created by pmat() function
#      with k = 3 and a fair coin.
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
        Pattern1 <- pattern_names[i],
        Pattern2 <- pattern_names[j],
        Expected_Wait = expected_wait[i],
        Win_Probability = M[i, j]
        )
      )
    }
  }
}

# Display pairs whose win probabilities notable differ from 0.5 
# (win probability is below 0.4 or above 0.6.)
print(comparison_table[
  abs(comparison_table$Win_Probability - 0.5) > 0.1, ])
  
  
  
