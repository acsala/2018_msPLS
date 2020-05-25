library(elasticnet)

calculateVectorEnet = function(x,y, lambda, nonzero){
  #function for elastic net penalization or multiple regression
  #
  #input:     x         (matrix - independent variables),
  #           y         (matrix - dependent variables),
  #           lambda    (integer - factor for Ridge penalty)
  #           nonzero   (integer - factor for LASSO penalty)
  
  #output:    tmp       vector    peanalized weights after regressing y on x
  
  nonzero = nonzero +1
  
  colnames(x) = paste("x", 1:ncol(x),sep=".")
  epsilon.enet = enet(x, y ,lambda=lambda, max.steps=nonzero)[[4]]
  
  #reorder
  tmp = rep(0,ncol(x))
  
  
  #Empty vector of regression coefficients
  tmp[colnames(x) %in% colnames(epsilon.enet)] = epsilon.enet[nonzero,]
  
  names(tmp) <- colnames(x)
  
  tmp
  
}