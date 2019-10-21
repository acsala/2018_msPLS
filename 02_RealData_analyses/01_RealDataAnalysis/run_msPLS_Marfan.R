load("/trinity/home/cscala/01attila/06_msPLS_Bioinformatics/02_RealData_analyses/01_RealDataAnalysis/datasetXYZ.RData")

# ALL 3 DATASETS####

data_sets <- cbind(X,Y,Z)

# path matrix
METHYL = c(0,1,0)
EXPRES = c(1,0,0)
CYTO = c(0,1,0)
sat.inner = rbind(METHYL, EXPRES, CYTO)

# blocks of outer model
sat.outer = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2], dim(Y)[2]+1:dim(Z)[2])
# define vector of reflective modes
sat.mod = c("B","B", "A")

library(sPLSPM)

time_data <- system.time(
  s_satpls2 <- splspm(data_sets, sat.inner, sat.outer, sat.mod, scheme="path",
                      scaled=T, penalization = "ust", 
                      nonzero = c(25,50), lambda = 1,
                      cross_validate = T)
)

res1.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "METHYL"),]
res2.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "EXPRES"),]
res3.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "CYTO"),]

res1.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="METHYL"),]
res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]
res3.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="CYTO"),]

# Sum abs weights of Y with latent variable of X
sum(abs(res2.inner[,"EXPRES"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls2$scores[,"EXPRES"],Y)))


# Sum abs weights of Z with latent variable of Y
sum(abs(res3.inner[,"EXPRES"]))
# Sum abs correlation of latent variable X and the Y variables
sum(abs(cor(s_satpls2$scores[,"EXPRES"],Z)))




# plot path coefficients
# plot(s_satpls, what="inner")

# plot loadings
# plot(s_satpls, what="loadings")
# 
# # plot outer weights
# plot(s_satpls, what="weights")

save(s_satpls2,time_data,cv_results, file = "sPLM_results.Rdata")