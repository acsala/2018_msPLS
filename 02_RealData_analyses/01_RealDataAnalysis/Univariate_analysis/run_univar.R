load("../datasetXYZ.RData")

library(xtable)


# FIRST DATASET
X[1:10, 1:10]

dim(X)[1]

#get first 5 patient 5 var
to_print <- X[1:5,1:5]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last two patients
to_print <- X[(dim(X)[1]-1):dim(X)[1],1:5]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last two vars
to_print <- X[1:5,(dim(X)[2]-1):dim(X)[2]]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last patient last var
to_print <- X[dim(X)[1],dim(X)[2]]
print(to_print, digits = 2)

# SECOND DATASET
#get first 5 patient 5 var
to_print <- Y[1:5,
              1:5]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last two patients
to_print <- Y[(dim(Y)[1]-1):dim(Y)[1],
              1:5]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last two vars
to_print <- Y[1:5,
              (dim(Y)[2]-1):dim(Y)[2]]
print(to_print, digits = 1)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last patient last var
to_print <- Y[dim(Y)[1],
              dim(Y)[2]]
print(to_print, digits = 2)


# THIRD DATASET
#get first 5 patient 5 var
to_print <- Z[1:5,
              1:5]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last two patients
to_print <- Z[(dim(Y)[1]-1):dim(Y)[1],
              1:5]
print(to_print, digits = 2)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last two vars
to_print <- Z[1:5,
              (dim(Z)[2]-1):dim(Z)[2]]
print(to_print, digits = 1)
print(xtable(to_print, type = "latex", digits = 2), include.rownames=FALSE)

#get last patient last var
to_print <- Z[dim(Z)[1],
              dim(Z)[2]]
print(to_print, digits = 2)


#*****************************************************************************************************
##Data analysis ####
#*****************************************************************************************************

# Pairwise

i <- 1
j <- 1

cor(X[,i],Y[,j])

X[,i] %*% Y[,j] / (dim(X)[1]-1)

(X[,i]-mean(X[,i])) %*% (Y[,j]-mean(Y[,j])) / 
  (( sqrt( sum( (X[,i]-mean(X[,i]) )^2) / (dim(X)[1]-1) )  ) %*%
     ( sqrt( sum( (Y[,j]-mean(Y[,j]) )^2) / (dim(Y)[1]-1) )  ) %*%
     (dim(X)[1]-1)
   )

lm(Y[,j]~X[,i])

plot(X[,i],Y[,j],
     xlab = colnames(X)[i],
     ylab = colnames(Y)[j])

abline(lm(Y[,j]~X[,i]) )

dev.off()

# Multivariable

i <- c(1,3)
j <- 1

#X[1:5,i]

var(X[,i]) 
#t(X[,i]) %*% X[,i] / (dim(X)[1]-1)

solve(var(X[,i])) %*% t(X[,i]) %*% Y[,j] / (dim(X)[1]-1)


lm(Y[,j]~X[,i])






to_print <- c()

to_print <- cbind(to_print,levels(names_char), 
                  rep("&",n), 
                  format((data_dealing_NOW$weight), digits = 2),
                  rep("&",n))

to_print <- cbind(to_print,levels(names_char2), 
                  rep("&",n), 
                  format((data_dealing_NOW2$weight), digits = 2),
                  rep("&",n))

to_print <- cbind(to_print,levels(names_char3), 
                  rep("&",n), 
                  format((data_dealing_NOW3$weight), digits = 2),
                  rep("\\",n))

to_print[,c(1,3)][6:30,] <- ""

print(data.frame(to_print),row.names = F)

X1 <- X[,1:36]

Betas <- solve(t(X1)%*%X1) %*% t(X1) %*% Y

absBetas <- apply(Betas,2,abs)

sum_abs <- apply(absBetas,2,sum)

sum_abs[order(sum_abs,decreasing=T)][1:100]
