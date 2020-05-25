##Peanlization algorithm with Univariate soft-thresholding - UST; ridge=Inf
get_ust = function(x,y, nonzero)
{
  beta = as.numeric(t(x)%*%y)
  
  #Select the highest values to be non-zero
  #*setting RIDGE to infinity indicates that highest values will be selected?
  return = rep(0, length(beta))
  
  select = order(abs(beta), decreasing=TRUE)[1:nonzero]
  
  return[select] = beta[select]
  
  names(return) <- colnames(x)
  
  return
}