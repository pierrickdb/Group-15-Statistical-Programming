#Add Names and Student ID here 
#Raneems Tasks:
#Fagans Tasks:
#Pierricks Tasks:
list.files()
setwd('C:/Users/s3015664/Downloads/Practical-1/Practical-1') ## comment out of submitted
music_raw <- readLines("palestrina.txt")

#Task 1 / preprocessing the Data
#Preprocess the music data 
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


#Task 2 / Create tokens and matrices
#Construct the vector of unique notes.
notes <- unique(music_clean)
#Find the index of each element in notes
tokens <- match(x = music_clean, table = notes)
#Create the function to build the matrix M
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

##############
M <- make_matrix(tokens = tokens, end = which(notes=="||"), mlag = 4)
head(M,100)
