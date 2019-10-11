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
library(xtable)
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

#univariable plot####

pdf("Univariate_plot_1.pdf") 

plot(X[,i],Y[,j],
     xlab = colnames(X)[i],
     ylab = colnames(Y)[j])

abline(lm(Y[,j]~X[,i]) )

dev.off()

# Multivariable

i <- 1:3
j <- 1

solve(var(X[,i])) %*% t(X[,i]) %*% Y[,j] / (dim(X)[1]-1)

# 
# solve(var(Y[,j])) %*%
#   solve(var(X[,i])) %*% t(X[,i]) %*% Y[,j] / (dim(X)[1]-1)

#X[1:5,i]

var(X[,i]) 
var(Y[,j])
#t(X[,i]) %*% X[,i] / (dim(X)[1]-1)
#print(xtable(var(X[,i]) , type = "latex", digits = 2), include.rownames=FALSE)

coefs <- solve(var(X[,i])) %*% t(X[,i]) %*% Y[,j] / (dim(X)[1]-1)
coefs

lm_model <- lm(Y[,j]~X[,i])

Y_hat <- X[,i] %*% coefs

Y[,j] - Y_hat

Y_hat
predict(lm_model)

#multivariable plot####

pdf("multivariable_plot_1.pdf") 
plot(X[,1],Y[,j],
     xlab = "Cpg sites",
     ylab = colnames(Y)[j])

for (vars in 1:dim(X[,i])[2]){
  
  points(X[,vars],Y[,j], col = vars)
  abline(lm(Y[,j] ~ X[,vars]), col = vars)
  
}

legend("topleft", legend=colnames(X)[i],
       col=c(1:3), 
       pch=1,
       #lty=1:2, 
       cex=0.8)
  
abline(Y_hat,Y[,j], lwd = 2)

dev.off()


#Run 2 dataset msPLS #####

data_sets <- cbind(X,Y)

# path matrix
METHYL = c(0,0)
EXPRES = c(1,0)
sat.inner = rbind(METHYL, EXPRES)

# blocks of outer model
sat.outer = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])
# define vector of reflective modes
sat.mod = c("B","B")

library(sPLSPM)

time_data <- system.time(
  s_satpls2 <- splspm(data_sets, sat.inner, sat.outer, sat.mod, scheme="path",
                      scaled=T, penalization = "ust", 
                      nonzero = c(37,37), lambda = 1,
                      cross_validate = F)
)

res1.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "METHYL"),]
res2.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "EXPRES"),]

res1.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="METHYL"),]
res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]

# Sum abs weights of Y with latent variable of X
sum(abs(res2.inner[,"EXPRES"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls2$scores[,"EXPRES"],Y)))


save(s_satpls2,time_data, file = "sCCA_res.RData.Rdata")


#plot 2 dataset pls lv results ####
Xi1 <- s_satpls2$scores[,1]
Xi2 <- s_satpls2$scores[,2]

pdf("CCA_LV_plot_1.pdf") 

plot(1, 
     type="n", 
     xlab=colnames(s_satpls2$scores)[1], 
     ylab=colnames(s_satpls2$scores)[2], 
     xlim=c(-3,3), 
     ylim=c(-3,3), 
     xaxt = "n",
     yaxt = "n",
     bty = "n")

points(Xi1,Xi2, 
       col = 1, 
       cex = 1, lwd = 1, pch = 1)


# legend("topleft", legend=colnames(X)[i],
#        col=c(1:3), 
#        pch=1,
#        #lty=1:2, 
#        cex=0.8)

abline(lm(Xi2~Xi1), 
       col = 1, lwd = 2)

axis(1, at = -3:3,
     labels = NULL,
     lwd = 2, 
     cex.axis = 1, las = 2)

axis(2, at = -3:3,
     labels = NULL,
     lwd = 2, 
     cex.axis = 1, las = 2)

#label points

pos_vector <- rep(3, length(Xi1))
pos_vector[names(Xi1) %in% c("183")] <- 1

text(Xi1, 
     Xi2, 
     #labels=names(Xi1), 
     labels=1:length(Xi1),
     cex= 1, pos=pos_vector)

dev.off()




nonzero_METHYL      <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="METHYL"),]
nonzero_EXPRESSION  <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="EXPRES"),]

plot_METHYL_alphas <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="METHYL"),]
plot_METHYL_alphas <- plot_METHYL_alphas[which(plot_METHYL_alphas[,3]!=0),][,c(1,3)]
plot_METHYL_alphas[]

plot_EXPRESSION_betas <- s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2] == "EXPRES"),]
plot_EXPRESSION_betas <- plot_EXPRESSION_betas[which(nonzero_EXPRESSION[,3]!=0),]
plot_EXPRESSION_betas









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
