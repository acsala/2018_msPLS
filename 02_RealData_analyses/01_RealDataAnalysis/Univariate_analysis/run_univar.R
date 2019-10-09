load("../datasetXYZ.RData")

X1 <- X[,1:36]

Betas <- solve(t(X1)%*%X1) %*% t(X1) %*% Y

absBetas <- apply(Betas,2,abs)

sum_abs <- apply(absBetas,2,sum)

sum_abs[order(sum_abs,decreasing=T)][1:100]
