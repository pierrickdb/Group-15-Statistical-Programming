#=================================================================================
#Add Raneems xy xx, Fagan xy xx, Pierrick de Boer s3015664
#Raneems Tasks: Task 1.4 Preprocessing
#Fagans Tasks: Task 5, 6 Matrix and Unique
#Pierricks Tasks: Task 1.1,1.2,1.3 7 & 8 Github, Read In, next_note, generate_music
#=================================================================================

#READ IN 
music_raw <- readLines("palestrina.txt")

#Input:
#Task:
#Return:
preprocess <- function(music) {
  #remove the first 3 lines
  music <- music[-c(1,2,3)] 
  #remove lines starting with #
  music <- music[-grep("^#", music)]
  #find the blank lines
  blank <- grep("^$", music)
  #add || to the end of each piece
  music[blank - 1] <- paste(music[blank - 1], "||")
  #add | to the end of the other measures
  music[-c(blank, blank - 1)] <- paste(music[-c(blank, blank - 1)],"|" )
  #remove the blank lines 
  music <- music[-blank]
  #split the music into individual notes
  music <- strsplit(music, " ") 
  #put everything into one vector 
  music <- unlist(music)
  return(music)
}

music_clean <- preprocess(music= music_raw)

#DICTIONARY AND TOKEN 
notes <- unique(music_clean)
tokens <- match(x = music_clean, table = notes)

#Input: List of tokens and token for end and size of considering window
#Task: 
#Return: Matrix M with all ...
make_matrix <- function(tokens, end, mlag = 4) {
  #Determine the dimension of matrix M
  n_column <- mlag + 1
  n_row <- length(tokens) - mlag
  #Make an empty vector for collecting data
  data <- c()
  #Fill vector with data from tokens
  for(i in 1:n_column) {
    index <- n_row - 1 + i
    data <- c(data, tokens[i:index]) 
  }
  #Create the matrix in column major order
  matrix <- matrix(data = data, nrow = n_row, ncol = n_column)
  #Remove all rows where end appears except the last column
  matrix <- matrix[-which(matrix[, 1:mlag] == end, arr.ind = TRUE)[, 1], ]
  return(matrix)
}

M <- make_matrix(tokens = tokens, end = which(notes=="||"), mlag = 4)

#Input: current key (vector of tokens regarding referring to the sequence), Matrix and weight
#Task: Creates the probability distribution which we will sample the next note from a given current sequence
#Return: the next_token
next_note <- function(key, M, w = rep(1,length(M[1,])-1)){
  
  #work with longer & shorter keys
  mlag <- length(M[1,]) - 1
  if (length(key) > mlag) {
    key <- key[(length(key) - mlag + 1):length(key)]
  }
  #start index
  mc <- mlag - length(key) + 1
  #storage for candidate and its probability
  cand = c()  
  prob = c()
  #consider all length form 1 to mlag keys for the next token
  while (mc <= mlag) {
    #
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == key)) 
    #store the next_token candidate where the key does not deviate
    u  <- M[ii == 0, mlag + 1]
    n  <- length(u)
    if (n > 0) {
      cand <- c(cand, u)
      prob <- c(prob, rep(w[mc] / n, n))
    }
    #shorten the key and increase the window
    key <- key[-1]   
    mc  <- mc + 1    
  }
  #sampling
  if (length(cand) == 0 || sum(prob) == 0) {
    next_token <- sample(M[,ncol(M)],1)
  } else if (length(cand) == 1) {
    next_token <- cand[1]                      
  } else {
    next_token <- sample(cand, 1, prob = prob)
  }
  return(next_token)
}

next_token <- next_note(key,M,w)

#Input: start token, M and weights
#Task: generate a vector with literately generating the next token until the song is complete
#Return the song as a concatenated string
generate_music <- function(start, M, w = rep(1, ncol(M) - 1) ,notes){
  key <- start
  amount_notes_used = 1
  amount_bars_used = 0
  #add notes until we have 200 notes or 8 bars or the song is ending (||)
  while(amount_notes_used < 200 && amount_bars_used <8){
    generated_token <- next_note(key,M,w)
    key <- c(key, generated_token)  
    amount_notes_used <- amount_notes_used + 1
    if(generated_token == which(notes == "|")) amount_bars_used <- amount_bars_used +1
    if(generated_token == which(notes == "||")) break 
    
  }
  #translate the tokens into actual notes and return song as concatenated string
  final_song <- notes[key]
  return(paste(final_song,collapse = " "))
}






