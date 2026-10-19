#=================================================================================
#Add Raneem Alluhaydan S3044929, Fagan Ganbarli S3015664, Pierrick de Boer S3034691
#Raneems Tasks: 
#Fagans Tasks: 
#Pierricks Tasks: Githubsetup Task 1 , 2, 3%
#=================================================================================


#install.packages("ggplot2") only if you did not install that yet 
library(ggplot2)

#x is the toss stream
#a is a pattern
firsthit <- function(x,a){
  
  #x must be a matrix now
  if(!is.matrix(x)){
    x <- matrix(x, nrow = 1)
  }
  
  L <- ncol(x)
  k <- length(a)
  
  #impossible case if the toss is smaller then pattern
  if(k>L){
    return(rep(NA, nrow(x)))
  }
  
  hit <- matrix(TRUE, nrow = nrow(x), ncol = L - k + 1)  
  for (j in seq_along(a)){
    hit <- hit & x[, j:(L - k + j), drop = FALSE] == a[j] 
  }
  

  #create return vector
  start <- max.col(hit, ties.method = "first")
  start[rowSums(hit) == 0] <- NA
  

  return(start + k - 1)
}

#returns a vector c(winne,tosses) where 1=player a and 0 player b and tosses is the number of duel 
#generate tosses in blocks, extending the stream only while neither pattern has yet appeared
sim_duels <- function(a,b,p=0.5,nreps = 1){
  
  k <- length(a)
  
  if(k!=length(b)){
    stop("this is not a fair game!")
  }
  if (all(a == b)) {
    stop("a and b are identical")
  }
  
  #intitial time and block length
  time <- 0
  B <- 3*k
  
  active <- 1:nreps
  winner <- rep(NA, nreps)
  tosses <- rep(NA, nreps)
  

  while(length(active) > 0){
    
    #create matrix block
    if(time == 0){
      W <- matrix(rbinom(nreps * B, size = 1, prob = p), nrow = nreps)
    }
    
    else{
      # letzte k-1 Spalten behalten (alle Zeilen)
      keep <- W[, B - (k - 1) + seq_len(k - 1), drop = FALSE]
      
      # neuer Block: für jedes aktive Duell B-k+1 neue Würfe
      W_new <- matrix(rbinom(nrow(keep) * (B - k + 1), size = 1, prob = p),
                      nrow = nrow(keep))
      
      # Spalten nebeneinander hängen
      W <- cbind(keep, W_new)
    }
    
    #check if pattern appears
    time_a <- firsthit(W,a) #(9,NA,10)
    time_b <- firsthit(W,b) #(1,8,12)
    
    a_wins <- !is.na(time_a) & (is.na(time_b) | time_a < time_b)#F,F,T
    b_wins <- !is.na(time_b) & (is.na(time_a) | time_b < time_a)#T,T,F
    
    winner[active[a_wins]] <- 1               
    winner[active[b_wins]] <- 0
    
    decided <- a_wins | b_wins
    hit_time <- pmin(time_a, time_b, na.rm = TRUE) 
    tosses[active[decided]] <- time + hit_time[decided]
    
    #now keep going with the rows that are not yes finished (NAvsNA)
    active <- active[!decided]
    #now we can drop the once which are terminated in W and start the loop again. 
    W <- W[!decided, , drop = FALSE] 

    
  
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
  
  


