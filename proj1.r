#Add Names and Student ID here 
#Raneems Tasks:
#Fagans Tasks:
#Pierricks Tasks:
list.files()
setwd('/Users/pierrickdeboer/Desktop/UoE-Master/Statistical Programming/Practical-1') ## comment out of submitted
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
