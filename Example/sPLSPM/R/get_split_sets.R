


get_split_sets <- function(X, Y, nr_subsets){
  #Function to split up data matrix into n subsets
  #
  #input:     X                   n*p matrix      - data matrix to be split
  #           Y                   n*q matrix      - data matrix to be split
  #           nr_subsets           integer         - nr of subsets the data matrix will be split to  
  #
  #
  #output:    label               vector with n_subset unique elements with length of n     
  #
  #           X.sampled,          n*p matrix      shuffled dataset of the origninal data
  #           Y.sampled,          n*q matrix      shuffled dataset of the origninal data
  
  
  n <- dim(X)[1]
  
  #re-sample data rows
  splitting_dimensions <- sample(1:n,n)
  
  X.sampled <- X[splitting_dimensions,]
  Y.sampled <- Y[splitting_dimensions,]
  
  #calculate how many rows are left over
  leftover = n %% nr_subsets
  
  rep_eat = (n-leftover)/nr_subsets 
  
  #repeat sequence nr_subsets 
  labels = rep(1:nr_subsets , each=rep_eat)
  
  if(leftover!=0)
  {
    labels = c(labels, (1:nr_subsets)[1:leftover])
  }
  
  #**********************
  result <- list(X.sampled = X.sampled,
                 Y.sampled = Y.sampled,
                 labels = labels
  )
  
  
  result
  
}