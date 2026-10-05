#Add Names and Student ID here 
#Raneems Tasks: Task 1.4 Preprocessing
#Fagans Tasks: Task 5, 6 Matrix and Unique
#Pierricks Tasks: Task 1.1,1.2,1.3 7 & 8 Github, readin, next_note, generate_music

#READ IN 
list.files()
setwd('C:/Users/s3015664/Downloads/Practical-1/Practical-1') ## comment out of submitted
music_raw <- readLines("palestrina.txt")

#TASK 1 PREPROCESSING
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
  # put everything into one vector 
  music <- unlist(music)
}
music_clean <- preprocess(music= music_raw)
tail(music_clean, n = 10)

#DICTIONARY AND TOKEN 
notes <- unique(music_clean)
tokens <- match(x = music_clean, table = notes)

#MATRIX 
make_matrix <- function(tokens, end, mlag = 4) {
  #Determine the dimension of matrix M
  n_column = mlag + 1
  n_row = length(tokens) - mlag
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
}

M <- make_matrix(tokens = tokens, end = which(notes=="||"), mlag = 4)
head(M,100)


#GENERATE NEXT NOTE
next_ <- function(key, M, w = rep(1,length(M[1,])-1)){
  
  #work with longer & shorter keys
  mlag <- length(M[1,]) - 1
  if (length(key) > mlag) {
    key <- key[(length(key) - mlag + 1):length(key)]
  }
  
  mc <- mlag - length(key) + 1
  cand = c()  #the next_note tokens stored
  prob = c()
  
  while (mc <= mlag) {
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == key)) #get the number of occurence
    u  <- M[ii == 0, mlag + 1]        # Treffer dieses Durchlaufs
    n  <- length(u)
    if (n > 0) {
      cand <- c(cand, u)
      prob <- c(prob, rep(w[mc] / n, n))
    }
    
    key <- key[-1]   # erstes Element weg
    mc  <- mc + 1    # Spaltenbereich wird entsprechend schmaler
  }
  
  #now we have both the prob and u and we can sample from them 
  if (length(cand) == 0 || sum(prob) == 0) {
    next_token <- sample(M[,ncol(M)],1)
  } else if (length(cand) == 1) {
    next_token <- cand[1]                      
  } else {
    next_token <- sample(cand, 1, prob = prob)
  }
  return(next_token)
}



