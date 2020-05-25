# Metadata ####
#
# Working example and simulation study for msPLS
#
# The R implementation of msPLS relies on the PLS-PM: Partial Least Squares Path Modeling R package from Sanchez et al. (2009).
# The package is modified so that it fits the requirements of big data analyis, as described in our manuscript.
# Below We provide an one of the data analysis studies we used in the manuscript. 

# Load the functions
library(R.utils)
library(mvtnorm)
library(tester)
library(turner)
library(diagram)
library(shape)
library(amap)
sourceDirectory("./sPLSPM/R/")

# replicatate the simulation study with n=250 ####

simulate_data <- function(sample_size = 100, 
                          nr_variables = c(50,50,20), 
                          nr_correlated_vars = c(5,5,5),
                          weight_correlations = c(0.7,0.5,0.3),
                          thetas = c(0.7,0.8,0.5),
                          sigma_correlation = 0.3){
  
  N = sample_size
  
  p1=nr_variables[1]
  p2=nr_variables[2]  
  p3=nr_variables[3]
  k=nr_correlated_vars[1]
  
  w1_cor = weight_correlations[1]
  w2_cor = weight_correlations[2]
  w3_cor = weight_correlations[3]
  
  theta1 = thetas[1]
  theta2 = thetas[2]
  
  Sigma = diag((p1+p2))
  Sigma[Sigma == 0] = 0
  #cor(X) 
  
  #Eerste k gecorreleerd met de eerste waarde uit Y
  Sigma[(1:k)+p1,(1:k)] = Sigma[(1:k),(1:k)+p1] = sigma_correlation
  
  data = rmvnorm(N, ,Sigma)
  
  X1 = data[,1:p1]
  X2 = data[,(p1+1):(p1+p2)]
  
  #corrplot(cor(X1[,1:10], X2[,1:10]))
  
  w1 = rep(0,p1)
  w1[1:k] = rep(w1_cor,k)
  
  w2 = rep(0,p2)
  w2[1:k] = rep(w2_cor,k)
  
  ksi1 = X1%*%w1
  ksi2 = X2%*%w2
  
  #solve(t(X1)%*%X1) %*% t(X1) %*% ksi1
  
  meanx = (theta1*ksi1) + (theta2*ksi2)
  sdx = (theta1)^2 + (theta2)^2
  
  ksi3 = rnorm(N,meanx,sqrt(abs(1-sdx)))
  
  w3 = rep(0,p3)
  w3[1:k] = rep(w3_cor,k)
  
  X3 = matrix(rep(0,N*p3),N,p3)
  
  for (j in 1:k) {
    sdx=sqrt(1-(w3)^2)
    X3[,j]=rnorm(N,meanx,sdx)      # sample values for the manifest variables
  }
  
  for (j in (k+1):p3) {
    X3[,j]=rnorm(N,0,1)      # sample values for the manifest variables
  }
  
  data$X1 <- X1
  data$X2 <- X2
  data$X3 <- X3
  
  
  return( data )
  
}

# Please change the sample size to replicate the other simulation study scenarios
sample_size = 250 
nr_variables = c(1000,1000,100)
nr_correlated_vars = c(10,10,10)
weight_correlations = c(0.7,0.6,0.3)
thetas = c(0.8,0.7,0.3)
sigma_correlation = 0.5
nonzero = c(5,10,15)

p0 = nr_variables - nr_correlated_vars
p1 = nr_correlated_vars
N = sample_size    # number of individuals
k = 3      # number of datasets
m = 1      # number of latent variables (LV's) per dataset

#repeat simulation ####

nr_of_simulations <- 5

sens.m <- matrix(c(0,0),nrow = max(nr_of_simulations),ncol = 2)
spec.m <- matrix(c(0,0),nrow = max(nr_of_simulations),ncol = 2)
ridge_param <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
lasso_param <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)

for (i in 1:nr_of_simulations){
  
  print("nr of simulation")
  print(i)    
  
  # for replication
  nr_seed = runif(1, 10, 10^8)
  set.seed(nr_seed)
  
  data <- simulate_data(sample_size = sample_size, 
                        nr_variables = nr_variables, 
                        nr_correlated_vars = nr_correlated_vars,
                        weight_correlations = weight_correlations,
                        thetas = thetas,
                        sigma_correlation = sigma_correlation)
  X1 <- data$X1
  X2 <- data$X2
  X3 <- data$X3
  
  Data <- cbind(X1,X2,X3)
  EXPL_X1 = c(0,1,0)
  EXPL_X2 = c(1,0,0)
  RESP_X3 = c(1,1,0)
  
  dim(Data)
  
  path_matrix = rbind(EXPL_X1, EXPL_X2, RESP_X3)
  #innerplot(path_matrix)
  
  # blocks of outer model
  blocks = list(1:dim(X1)[2], 
                (dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),
                (dim(X1)[2]+dim(X2)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]))
  
  modes = c("B","B","A")
  
  time_data <- system.time(
    s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                       scaled=T, penalization = "ust", nonzero = nonzero, 
                       lambda = 1, maxiter = 100, cross_validate = T)
  )
  
  print("s_satpls$model$iter")
  print(s_satpls$model$iter)
  
  print("s_satpls$nonzero")
  print(s_satpls$nonzero)
  print("s_satpls$lambda")
  print(s_satpls$lambda)
  
  ridge_param[i] <- s_satpls$lambda
  lasso_param[i] <- s_satpls$nonzero
  
  
  #calculate sens spec
  nzero_X_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X1",3])>0)
  nzero_Y_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X2",3])>0)
  
  zero_X_positiong <- which((s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X1",3])==0)
  zero_Y_positiong <- which((s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X2",3])==0)
  
  #X sensitivity
  sensX = sum(nzero_X_positiong %in% 1:nr_correlated_vars[1])/min(nr_correlated_vars[1],length(nzero_X_positiong))
  sensY = sum(nzero_Y_positiong %in% 1:nr_correlated_vars[2])/min(nr_correlated_vars[2],length(nzero_Y_positiong))
  
  print("print(c(sensX,sensY))")
  print(c(sensX,sensY))
  
  
  #specificiy
  specX = sum(zero_X_positiong %in% (nr_correlated_vars[1]+1):nr_variables[1])/(nr_variables[1]-nr_correlated_vars[1])
  specY = sum(zero_Y_positiong %in% (nr_correlated_vars[2]+1):nr_variables[2])/(nr_variables[2]-nr_correlated_vars[2])
  
  print("print(c(specX,specY))")
  print(c(specX,specY))
  
  
  sens.m[i,]  =c(sensX,sensY)
  spec.m[i,]  =c(specX,specY)
  
  
  iter    = s_satpls$model$iter
  nonzero = s_satpls$nonzero
  lambda  = s_satpls$lambda
  
  
  print("matrix Sens:")
  print(sens.m)
  print(apply(sens.m,2,mean))
  print("matrix spec:")
  print(spec.m)
  print(apply(spec.m,2,mean))
  
}

print("Mean Sens:")
print(apply(sens.m,2,mean))
print("Mean spec:")
print(apply(spec.m,2,mean))