case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jeremyrcoyle/gentmle2")
# devtools::install_github("jlstiles/Simulations", force = TRUE)
library(Simulations)
source("WrappersVblip1.R")

n = 2000
B = 100

dgps = lapply(1:B, FUN = function(x) get.dgp(n,4))

detectCores()
cl = makeCluster(12, type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)


# debug(SL.stack1)
# debug(sim_cv)
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.hal")
# SL.library = list("SL.nnet", "glm.mainint", c("SL.hal", "screen.Main"))
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.glm.interaction", "SL.mean","SL.rpartPrune", "SL.earth", "SL.glmnet")
# SL.library = list("SL.nnet", "glm.mainint", "SL.mean", "SL.glm","SL.rpartPrune", "earthFull","SL.glmnet")
# SL.library = SL.libraryG = c("SL.glm", "SL.mean")

simHal = function(data, gform = gform, Qform = Qform, 
                  V = 10, single = FALSE, estimator, method = "method.NNloglik", 
                  gn = NULL, family = binomial(), dgp) {
  
  S = sim_hal(data, gform = gform, Qform = Qform, 
          V = 10, single = FALSE, estimator, method, 
          gn = NULL, family, cvhal = TRUE)
  res = c(S, BV0 = dgp$BV0, ATE0 = dgp$ATE0)
  
  return(list(res = res, gn = dgp$PGn, blip_n = dgp$blip_n))
}
gform = formula("A~.")
Qform = formula("Y~A*(W1+W2+W3+W4)")
ALL=foreach(i=1:1,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {simHal(data = dgps[[i]]$DF, gform = gform, Qform = Qform, 
                    V = 10, single = FALSE, estimator = "single 1step", method = "method.NNloglik", 
                    gn = NULL, family = binomial(), dgp = dgps[[i]]) 
              }
ALL
save(ALL, dgps, file = "caseHalvsDelta.RData")

# res = lapply(ALL, FUN = function(x) x$res)
# results = data.matrix(data.frame(do.call(rbind, res)))


# colnames(results)
# 
# ateInds
# bvInds
# 
# vapply(ateInds, FUN = function(x) {
#   mean(results[,x+1] <= results[,"ATE0"] & results[,x+2] >= results[,"ATE0"])
# }, FUN.VALUE = 1)
# 
# vapply(bvInds, FUN = function(x) {
#   mean(results[,x+1] <= results[,"BV0"] & 
#          results[,x+2] >= results[,"BV0"])
# }, FUN.VALUE = 1)
# 
# bvAll
# ateAll
# 
# vapply(bvAll, FUN = function(x) {
#   mean(results[,x] - results[, "BV0"])
# }, FUN.VALUE = 1)
# 
# vapply(ateAll, FUN = function(x) {
#   mean(results[,x] - results[,"ATE0"])
# }, FUN.VALUE = 1)
# 
# BIAS = vapply(ateAll, FUN = function(x) {
#   results[,x] - results[,"ATE0"]
# }, FUN.VALUE = rep(1, nrow(results)))
# 
# hist(BIAS[,1])
# 
# vapply(ateInds, FUN = function(x) {
#   mean(results[,x+2] - results[,x+1])
# }, FUN.VALUE = 1)
# 
# vapply(bvInds, FUN = function(x) {
#   mean(results[,x+2] - results[,x+1])
# }, FUN.VALUE = 1)
# #
# colMeans(results[,grep("coef", colnames(results))])

#
#
# results[,c(4, 24, 30, 42)]
# for (i in 1:nrow(results)) hist(dgps[[i]]$blip_n, breaks = 200)
#
# for (i in 1:nrow(results)) hist(dgps[[i]]$PQ1n, breaks = 200)
# for (i in 1:nrow(results)) hist(dgps[[i]]$PQ0n, breaks = 200)
#
# vapply(1:nrow(results), FUN = function(x) {
#   c(dgps[[x]]$BV0, dgps[[x]]$ATE0)
# }, FUN.VALUE = c(1,1))
#
# vapply(1:nrow(results), FUN = function(x) {
#   c(min(dgps[[x]]$PGn), max(dgps[[x]]$PGn))
# }, FUN.VALUE = c(1,1))
#

# 
# i = 4
# YA1 = mean((dgps[[i]]$DF$Y/dgps[[i]]$PGn)[dgps[[i]]$DF$A==1])
# YA0 = mean((dgps[[i]]$DF$Y/(1-dgps[[i]]$PGn))[dgps[[i]]$DF$A==0])
# YA1
# YA0
# YA1 - YA0
# dgps[[i]]$ATE0
# mean(dgps[[i]]$blip_n)
