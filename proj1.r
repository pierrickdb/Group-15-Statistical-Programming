#=================================================================================
#Add Raneems S3044929, Fagan S3015664, Pierrick de Boer S3034691
#Raneems Tasks: Task 1.4 Preprocessing
#Fagans Tasks: Task 5, 6 Matrix and Unique
#Pierricks Tasks: Task 1.1,1.2,1.3 7 & 8 Github, Read In, next_note, generate_music
#=================================================================================

#READ IN 
music_raw <- readLines("palestrina.txt")

#Input: music data
#Task: preprocess the music data
#Return: cleaned music data

preprocess <- function(music) {
  
  #remove the first 3 lines
  music <- music[-c(1,2,3)] 
  
  #remove lines starting with #
  music <- music[-grep("^#", music)]
  
  #find the number of lines
  n <- length(music)
  
  # go through each line 
  for(i in 1:(n-1)) {
    
    #add || before a blank line
    
    if (music[i+1] == "") {  
      music[i] <- paste(music[i], "||")
      
      
      
    }else{
      #add | to the other lines
      music[i] <- paste(music[i], "|")
    }
  }
  
  
  #remove blank lines 
  music <- music[music != ""]
  
  #split the music into individual notes
  music <- strsplit(music, " ") 
  
  # put everything into one vector 
  music <- unlist(music)
  
  return(music)
}

music_clean <- preprocess(music= music_raw)

tail(music_clean, n = 10)


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
  remove_rows <- which(rowSums(M[, 1:mlag] == end) > 0)
  if (length(remove_rows) > 0) {
    M <- M[-remove_rows, , drop = FALSE]
  }
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

start <- 20
w3 <- c(0, 0, 0, 1)
final_music <- generate_music(start,M,w3,notes)
final_music

# Experiment with different mixture weights
#
# w1 = c(1,1,1,1):
# All four lagged dependencies are used. The model therefore considers
# both long-term and short-term patterns in the music. This provides
# more context when predicting the next note and can help preserve
# stylistic features such as note sequences and bar structure.
#
# w2 = c(1,0,0,0):
# Only the longest lagged dependency is used. The prediction is based
# on the previous four notes, corresponding to a fourth-order Markov
# model. This gives the model more context, but matching four-note
# sequences occurs less frequently, so the predictions can be more
# variable when there are few matching sequences.
#
# w3 = c(0,0,0,1):
# Only the shortest lagged dependency is used. The prediction therefore
# depends only on the immediately preceding note, corresponding to a
# first-order Markov model. This gives less contextual information and
# therefore tends to produce less structured sequences, although local
# note transitions are still preserved.







