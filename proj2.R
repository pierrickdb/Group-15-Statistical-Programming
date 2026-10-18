#=================================================================================
#Add Raneem Alluhaydan S3044929, Fagan Ganbarli S3015664, Pierrick de Boer S3034691
#Raneems Tasks: 
#Fagans Tasks: 
#Pierricks Tasks: Githubsetup xy%
#=================================================================================


#install.packages("ggplot2") only if you did not install that yet 
library(ggplot2)

#x is the toss stream
#a is a pattern
firsthit <- function(x,a){
  
  #impossible case if the toss is smaller then pattern
  if(length(a)>length(x)){
    return(NA)
  }
  
  for (i in seq_along(x)){
      #find start
      if(x[i]==a[1]){
        #now all values have to be the same
        
        #edge case
        if (i + length(a) - 1 > length(x)) return(NA)
        
        for (j in seq_along(a)){

          #same we look at the next
          if(x[i]==a[j] && j == length(a)){
            return(i)
          }
          #all values are the same
          else if(x[i]==a[j]){
            i <- i + 1
          }
          
          #if not search the next strating value of that patter
          else{
            break
          }
          
        }
        
      }
  }

  #not found 
  return(NA)
}

#returns a vector c(winne,tosses) where 1=player a and 0 player b and tosses is the number of duel 
#generate tosses in blocks, extending the stream only while neither pattern has yet appeared
duel <- function(a,b,p=0.5){
  
  if(length(a)!=length(b)){
    stop("this is not a fair game!")
  }
  
  #intitial time and block length
  time = 0
  B<- 3*length(a)
  
  while(TRUE){
    
    #create block
    if(time == 0){
      block <- rbinom(n = B, size = 1, prob = p)
    }
    else{
      block <- tail(block,length(a)-1)
      block_new <- rbinom(n = B - length(a) + 1, size = 1, prob = p)
      block <- c(block,block_new)
      time <- time - (length(a)-1)
    }
    
    #check if pattern appears
    time_a <- firsthit(block,a)
    time_b <- firsthit(block,b)
    if((!is.na(time_a) && is.na(time_b)) || (!is.na(time_a) && !is.na(time_b) && time_a < time_b)){
      time_winner <- time_a + time
      winner <- 1
      return(c(winner,time_winner))
    }
    else if((is.na(time_a) && !is.na(time_b)) || (!is.na(time_a) && !is.na(time_b) && time_a > time_b)){
      time_winner <- time_b + time
      winner <- 0
      return(c(winner,time_winner))
    }
    else{
      #new block and add block length to the time
      time <- time + B
    }
  }

}


