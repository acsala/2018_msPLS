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

plotting_function <- function(lim_x = c(-3,3),
                              lim_y = c(-3,3),
                              x_name = "x axis",
                              y_name = "y axis",
                              main_name = "main"){
  
  
  plot(1, 
       main = main_name,
       type="n", 
       xlab=paste0(x_name), 
       ylab=paste0(y_name), 
       xlim=lim_x, 
       ylim=lim_y, 
       xaxt = "n",
       yaxt = "n",
       bty = "n")
  
  
  # legend("topleft", legend=colnames(X)[i],
  #        col=c(1:3), 
  #        pch=1,
  #        #lty=1:2, 
  #        cex=0.8)
  
  
  axis(1, at = lim_x,
       labels = NULL,
       lwd = 2, 
       cex.axis = 1, las = 2)
  
  axis(2, at = lim_y,
       labels = NULL,
       lwd = 2, 
       cex.axis = 1, las = 2)
  
}


plotting_function(main_name = bquote(xi[.(j)]),
                  x_name = paste0("Methylation site: ",colnames(X)[i]),
                  y_name = paste0("Gene expression: ",colnames(Y)[j]))

points(X[,i],Y[,j],
       col = 1, 
       cex = 1, lwd = 1, pch = 1)


abline(lm(Y[,j]~X[,i]), 
       col = "red",  
       lwd=3)


text(X[,i],Y[,j],
     #labels=names(Xi1), 
     labels=1:length(X[,i]),
     cex= 1, pos=4)

pdf("Univariate_plot_1.pdf") 

plot(X[,i],Y[,j], 
     main = bquote(hat(y[1]) == .(format(lm(Y[,j]~X[,i])$coef[2], digits=2)) ~ x[1] ),
     col.main = "red",
     xlab = paste0("Methylation site: ",colnames(X)[i]),
     ylab = paste0("Gene expression: ",colnames(Y)[j]))

abline(lm(Y[,j]~X[,i]), col = "red",  lwd=3)

text(X[,i],Y[,j],
     #labels=names(Xi1), 
     labels=1:length(X[,i]),
     cex= 1, pos=4)



format(lm(Y[,j]~X[,i])$coef[2], digits=2)

dev.off()


pdf("Univariate_plot_2.pdf") 

par(mfrow=c(2,2))

for(which_x in 1:4){
  
  plot(X[,which_x],Y[,j],
       main = bquote(hat(y[1]) == .(format(lm(Y[,j]~X[,which_x])$coef[2], digits=2)) ~ x[.(which_x)] ),
       col.main = "red",
       xlab = colnames(X)[which_x],
       ylab = colnames(Y)[j])
  
  abline(lm(Y[,j]~X[,which_x]), col = "red",  lwd=3)
  
}

dev.off()

par(mfrow=c(1,1))

## GGPlot example with labels ####
pdf("try.pdf") 
library(ggrepel)
set.seed(42)

dat <- subset(mtcars, wt > 2.75 & wt < 3.45)
dat$car <- rownames(dat)

p <- ggplot(dat, aes(wt, mpg, label = car)) +
  geom_point(color = "red")

p1 <- p + geom_text() + labs(title = "geom_text()")

p2 <- p + geom_text_repel() + labs(title = "geom_text_repel()") + theme_classic()

gridExtra::grid.arrange(p1, p2, ncol = 2)


set.seed(42)

dat2 <- subset(mtcars, wt > 3 & wt < 4)
# Hide all of the text labels.
dat2$car <- ""
# Let's just label these items.
ix_label <- c(2,3,16)
dat2$car[ix_label] <- rownames(dat2)[ix_label]

ggplot(dat2, aes(wt, mpg, label = car)) +
  geom_point(color = ifelse(dat2$car == "", "grey50", "red")) +
  geom_text_repel()

dev.off()
##

# Multivariable ######

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
summary(lm_model)

Y_hat <- X[,i] %*% coefs

Y_hat
predict(lm_model)

X[1,i] %*% coefs

#multivariable plot####

pdf("multivariable_plot_1.pdf") 
plot(X[,1],Y[,j],
     xlab = "Cpg sites",
     ylab = colnames(Y)[j])

for (vars in 1:dim(X[,i])[2]){
  
  points(X[,vars],Y[,j], col = vars)
  abline(lm(Y[,j] ~ X[,vars]), col = vars)
  
  # text(X[,vars],Y[,j],
  #      #labels=names(Xi1), 
  #      labels=names(X[,vars]),
  #      cex= 1, pos=4)
  
}

legend("topleft", legend=colnames(X)[i],
       col=c(1:3), 
       pch=1,
       #lty=1:2, 
       cex=0.8)
  
abline(Y_hat,Y[,j], lwd = 2)

dev.off()


pdf("multivariable_plot_2.pdf") 
plot(Y_hat,Y[,j],
     main = bquote(hat(y)[1] == .(format(lm(Y[,j]~X[,i])$coef[2], digits=2)) ~ x[1] ~ 
                     .(format(lm(Y[,j]~X[,i])$coef[3], digits=2)) ~ x[2] ~ 
                     .(format(lm(Y[,j]~X[,i])$coef[4], digits=2)) ~ x[3]  ),
     col.main = "red",
     xlim = c(-2,2),
     ylim = c(-2,2),
     xlab = "Linear combination of Cpg sites",
     ylab = colnames(Y)[j])

text(Y_hat,Y[,j],
     #labels=names(Xi1),
     labels=1:length((Y[,j])),
     cex= 1, pos=4)

abline(Y_hat,Y[,j], lwd = 2, col = "red")

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
                      nonzero = c(10,10), lambda = 1,
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
     xlab=paste0("Linear combination of the methylaton sites"), 
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



# print data we found with CCA #####
library(xtable)
nonzero_METHYL      <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="METHYL"),]
nonzero_EXPRESSION  <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="EXPRES"),]

plot_METHYL_alphas <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="METHYL"),]
plot_METHYL_alphas <- plot_METHYL_alphas[which(plot_METHYL_alphas[,3]!=0),][,c(1,3)]
plot_METHYL_alphas[]

plot_METHYL_betas <- s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2] == "METHYL"),]
plot_METHYL_betas <- plot_METHY_betas[which(nonzero_METHYL[,3]!=0),]
plot_METHYL_betas


plot_EXPRESSION_alphas <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2]=="EXPRES"),]
plot_EXPRESSION_alphas <- plot_EXPRESSION_alphas[which(plot_EXPRESSION_alphas[,3]!=0),][,c(1,3)]
plot_EXPRESSION_alphas

plot_EXPRESSION_betas <- s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2] == "EXPRES"),]
plot_EXPRESSION_betas <- plot_EXPRESSION_betas[which(nonzero_EXPRESSION[,3]!=0),]
plot_EXPRESSION_betas

cbind(as.character(plot_EXPRESSION_betas[,"name"]), format(plot_METHYL_betas[,"EXPRES"], digits = 2),
as.character(plot_METHY_betas[,"name"]), format(plot_EXPRESSION_betas[, "METHYL"], digits = 2))

print(
  xtable(cbind(as.character(plot_METHYL_betas[,"name"]), format(plot_METHYL_betas[,"EXPRES"], digits = 2),
               as.character(plot_EXPRESSION_betas[,"name"]), format(plot_EXPRESSION_betas[, "METHYL"], digits = 2)), 
         type = "latex", digits = 2), 
  include.rownames=FALSE)





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

# take second LV set ####
load("../datasetXYZ.RData")
load("sCCA_res.RData")

zeta1 <- s_satpls2$scores[,1]
zeta2 <- s_satpls2$scores[,2]

get_residuals <- function(Dataset, LV){
  
  # calculate the residuals
  calcres = function(Xcol)
    Xcol - solve(t(LV)%*%LV) %*% t(LV) %*% Xcol %*% t(LV)
  
  Res_data = apply(Dataset, 2, calcres)
  
  return(Res_data)
  
}

Res_X <- get_residuals(X, zeta1)
Res_Y <- get_residuals(Y, zeta2)

data_sets <- cbind(Res_X,Res_Y)

print(dim(data_sets))

# path matrix
METHYL = c(0,0)
EXPRES = c(1,0)
sat.inner = rbind(METHYL, EXPRES)

# blocks of outer model
sat.outer = list(1:dim(Res_X)[2], dim(Res_X)[2]+1:dim(Res_Y)[2])
# define vector of reflective modes
sat.mod = c("B","B")

library(sPLSPM)

time_data <- system.time(
  s_satpls2_res <- splspm(data_sets, sat.inner, sat.outer, sat.mod, scheme="path",
                      scaled=T, penalization = "ust", 
                      nonzero = c(10,10), lambda = 1,
                      cross_validate = F)
)

res1.outer_residual <- s_satpls2_res$outer_model[
  which(s_satpls2_res$outer_model[,2] == "METHYL"),]
res2.outer_residual <- s_satpls2_res$outer_model[
  which(s_satpls2_res$outer_model[,2] == "EXPRES"),]

res1.inner_residual <-  s_satpls2_res$crossloadings[
  which(s_satpls2_res$crossloadings[,2]=="METHYL"),]
res2.inner_residual <-  s_satpls2_res$crossloadings[
  which(s_satpls2_res$crossloadings[,2]=="EXPRES"),]

# Sum abs weights of Y with latent variable of X
sum(abs(res2.inner_residual[,"EXPRES"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls2_res$scores[,"EXPRES"],Y)))


save(s_satpls2_res,time_data, file = "sCCA_secondLVs_results.RData")

#plot 1st and 2nd components ####
load("sCCA_res.RData")
load("sCCA_secondLVs_results.RData")

zeta1 <- s_satpls2$scores[,1]
zeta2 <- s_satpls2$scores[,2]

zeta1_residual_1 <- s_satpls2_res$scores[,1]
zeta2_residual_1 <- s_satpls2_res$scores[,2]

cor(zeta1, zeta1_residual_1)


pdf("CCA_LV1_vsLV2_plot_1.pdf") 

plot(1, 
     type="n", 
     xlab=paste0("First component scores (13%)"), 
     ylab=paste0("Second component scores (61%)"), 
     xlim=c(-3,3), 
     ylim=c(-3,3), 
     xaxt = "n",
     yaxt = "n",
     bty = "n")

points(zeta1,zeta1_residual_1, 
       col = 1, 
       cex = 1, lwd = 1, pch = 1)


# legend("topleft", legend=colnames(X)[i],
#        col=c(1:3), 
#        pch=1,
#        #lty=1:2, 
#        cex=0.8)

# abline(lm(zeta1_residual_1~zeta1), 
#        col = 1, lwd = 2)

axis(1, at = -3:3,
     labels = NULL,
     lwd = 2, 
     cex.axis = 1, las = 2)

axis(2, at = -3:3,
     labels = NULL,
     lwd = 2, 
     cex.axis = 1, las = 2)

#label points

pos_vector <- rep(3, length(zeta1))
pos_vector[names(zeta1) %in% c("183")] <- 1

text(zeta1, 
     zeta1_residual_1, 
     #labels=names(Xi1), 
     labels=1:length(zeta1),
     cex= 1, pos=pos_vector)

dev.off()

#plot first and second component loading dimensions ####

plot_METHYL <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "METHYL"),]
plot_METHYL <- plot_METHYL[which(plot_METHYL[,"weight"]!=0),][,]
plot_METHYL$loading

plot_METHYL_resid1 <- s_satpls2_res$outer_model[
  which(s_satpls2_res$outer_model[,2] == "METHYL"),]
plot_METHYL_resid1 <- plot_METHYL_resid1[which(plot_METHYL_resid1[,"weight"]!=0),][,]
plot_METHYL_resid1$loading

plotting_LVs_n_scores <- function(X_plot, Y_plot,
                                  X_lim, Y_lim,
                                  X_text, Y_text){
  
  X_lim = c(min(X_lim)-sd(X_lim),
            max(X_lim)+sd(X_lim))
  Y_lim = c(min(Y_lim)-sd(Y_lim),
            max(Y_lim)+sd(Y_lim))

  plot(1, 
       type="n", 
       xlab=X_text, 
       ylab=Y_text, 
       xlim=X_lim, 
       ylim=Y_lim, 
       xaxt = "n",
       yaxt = "n",
       bty = "n",
       lwd = 2)
  
  points(X_plot,Y_plot, 
         col = 1, 
         cex = 1, lwd = 1, pch = 1)
  
  
  # legend("topleft", legend=colnames(X)[i],
  #        col=c(1:3), 
  #        pch=1,
  #        #lty=1:2, 
  #        cex=0.8)
  
  # abline(lm(zeta1_residual_1~zeta1), 
  #        col = 1, lwd = 2)
  
  axis(1, 
       at = format(seq(from = X_lim[1], to=X_lim[2], by = (X_lim[2]-X_lim[1])/5), 
                      digits=2),
       labels = NULL,
       lwd = 2, 
       cex.axis = 1, las = 2)
  
  axis(2, at = format(seq(from = Y_lim[1], to=Y_lim[2], by = (Y_lim[2]-Y_lim[1])/5), 
                      digits=2),
       labels = NULL,
       lwd = 2, 
       cex.axis = 1, las = 2)
  
}

Methyl_weights_cor_1st <-
  cor(X[,colnames(X) %in% plot_METHYL$name | 
          colnames(X) %in% plot_METHYL_resid1$name],
      s_satpls2$scores[,1])

Methyl_weights_cor_2nd <-
  cor(X[,colnames(X) %in% plot_METHYL$name | 
          colnames(X) %in% plot_METHYL_resid1$name],
      s_satpls2_res$scores[,1])

cor(Methyl_weights_cor_1st, Methyl_weights_cor_2nd)

pdf("CCA_METHYL_1st_2nd_component.pdf") 

plotting_LVs_n_scores(Methyl_weights_cor_1st, Methyl_weights_cor_2nd,
                      X_lim = range(Methyl_weights_cor_1st),
                      Y_lim = range(Methyl_weights_cor_2nd),
                      X_text = "1st component loadings on CpG sites",
                      Y_text = "2nd component loadings on CpG sites"
                      )


text(Methyl_weights_cor_1st[
  rownames(Methyl_weights_cor_1st) %in% plot_METHYL_resid1$name],
     1.5,
     #labels=names(Xi1),
     labels=rownames(Methyl_weights_cor_1st)[
         rownames(Methyl_weights_cor_1st) %in% plot_METHYL_resid1$name],
     cex= 1,
     srt = 90, 
  col = 1)

text(1.2,
     Methyl_weights_cor_2nd[
       rownames(Methyl_weights_cor_2nd) %in% plot_METHYL$name],
  #labels=names(Xi1),
  labels=rownames(Methyl_weights_cor_2nd)[
    rownames(Methyl_weights_cor_2nd) %in% plot_METHYL$name],
  cex= 1,
  srt = 0, col = 1)

dev.off()

plot_EXPRESSION <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "EXPRES"),]
plot_EXPRESSION <- plot_EXPRESSION[which(plot_EXPRESSION[,3]!=0),][,]
plot_EXPRESSION$loading

plot_EXPRESSION_resid1 <- s_satpls2_res$outer_model[
  which(s_satpls2_res$outer_model[,2] == "EXPRES"),]
plot_EXPRESSION_resid1 <- plot_EXPRESSION_resid1[which(plot_EXPRESSION_resid1[,3]!=0),][,]
plot_EXPRESSION_resid1$loading

Expr_weights_cor_1st <-
  cor(Y[,colnames(Y) %in% plot_EXPRESSION$name | 
          colnames(Y) %in% plot_EXPRESSION_resid1$name],
      s_satpls2$scores[,1])

Expr_weights_cor_2nd <-
  cor(Y[,colnames(Y) %in% plot_EXPRESSION$name | 
          colnames(Y) %in% plot_EXPRESSION_resid1$name],
      s_satpls2_res$scores[,1])


pdf("CCA_EXPR_1st_2nd_component.pdf") 

plotting_LVs_n_scores(Expr_weights_cor_1st, Expr_weights_cor_2nd,
                      X_lim = range(Expr_weights_cor_1st),
                      Y_lim = range(Expr_weights_cor_2nd),
                      X_text = "1st component loadings on Expressions",
                      Y_text = "2nd component loadings on Expressions"
)


text(Expr_weights_cor_1st[
  rownames(Expr_weights_cor_1st) %in% plot_EXPRESSION_resid1$name],
  -1,
  #labels=names(Xi1),
  labels=rownames(Expr_weights_cor_1st)[
    rownames(Expr_weights_cor_1st) %in% plot_EXPRESSION_resid1$name],
  cex= 1,
  srt = 90, 
  col = "red")

points(Expr_weights_cor_1st[
  rownames(Expr_weights_cor_1st) %in% plot_EXPRESSION_resid1$name],
  Expr_weights_cor_2nd[
    rownames(Expr_weights_cor_1st) %in% plot_EXPRESSION_resid1$name], 
       col = "red", 
       cex = 1, lwd = 1, pch = 1)

text(0,
     Expr_weights_cor_2nd[
       rownames(Expr_weights_cor_2nd) %in% plot_EXPRESSION$name],
     #labels=names(Xi1),
     labels=rownames(Expr_weights_cor_2nd)[
       rownames(Expr_weights_cor_2nd) %in% plot_EXPRESSION$name],
     cex= 1,
     srt = 0, col = "blue")

points(Expr_weights_cor_1st[
  rownames(Expr_weights_cor_2nd) %in% plot_EXPRESSION$name],
  Expr_weights_cor_2nd[
    rownames(Expr_weights_cor_2nd) %in% plot_EXPRESSION$name], 
  col = "blue", 
  cex = 1, lwd = 1, pch = 1)

dev.off()
