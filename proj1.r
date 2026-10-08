#=================================================================================
#Add Raneems S3044929, Fagan S3015664, Pierrick de Boer S3034691
#Raneems Tasks: Task 1.4 Preprocessing
#Fagans Tasks: Task 5, 6 Matrix and Unique
#Pierricks Tasks: Task 1.1,1.2,1.3 7 & 8 Github, Read In, next_note, generate_music
#=================================================================================

#READ IN
music_raw <- readLines("palestrina.txt")

#Input: raw music data
#Task: preprocess the music data
#Return: cleaned music data
preprocess <- function(music) {
  
  #remove the first 3 lines
  music <- music[-c(1,2,3)] 
  
  #remove lines starting with #
  music <- music[-grep("^#", music)]
  
  #find the number of lines
  n <- length(music)
  
  #go through each line 
  for(i in 1:(n-1)) {


    if (music[i] == "") {
      #skip blank lines to avoid adding | to their end
      next
      
    }else if (music[i+1] == "") {  
      #add || before a blank line
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

#Input: preprocessed music 
#Task: create a vector of unique notes
#Return: the notes(dictionary)
notes <- unique(music_clean)

#Input: preprocessed music, unique notes
#Task: create a vector of token IDs that refers to each element of music_clean
#Return: list of tokens
tokens <- match(x = music_clean, table = notes)

#Input: list of tokens, mlag which is the amount of previous notes for determining the next note 
#Task: create the function to build the matrix M
#Return: Matrix M containing token sequences
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
  matrix_ <- matrix(data = data, nrow = n_row, ncol = n_column)
  
  #Remove rows where "end" appears in the first mlag columns
  matrix_ <- matrix_[-which(matrix_[, 1:mlag] == end, arr.ind = TRUE)[, 1], ]

  return(matrix_)
}

M <- make_matrix(tokens = tokens, end = which(notes=="||"), mlag = 4)

#Input: current key,  Matrix and weight
#Task: creates the probability distribution which we will sample the next note from a given current sequence
#Return: the following token 
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
    
    #reapeatedly finding the rows of M that matches the key 
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


#Input: start token, M and weights
#Task: generate a vector with iteratively calling the next token and decode it into notes
#Return the song as a concatenated string
generate_music <- function(start, M, w = rep(1, ncol(M) - 1) ,notes){

  key <- start
  amount_notes_used = 1
  amount_bars_used = 0
  
  #add notes until we have 200 notes, 8 bars or the song is ending (||)
  while(amount_notes_used < 200 && amount_bars_used <8){
    
    generated_token <- next_note(key,M,w)
    
    key <- c(key, generated_token)  
    
    amount_notes_used <- amount_notes_used + 1
    
    #count bars
    if(generated_token == which(notes == "|")) {
      amount_bars_used <- amount_bars_used +1
    }  
    #stop sign
    if(generated_token == which(notes == "||")) {
      break 
    }
    
  }
  #translate token into note 
  final_song <- notes[key]
  
  return(paste(final_song,collapse = " "))
}


sim_music <- generate_music(start = 20, M = M, notes = notes,
                            w = c(1, 1, 1, 1))
#print song
sim_music 


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







