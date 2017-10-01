case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jeremyrcoyle/gentmle2")
# devtools::install_github("jlstiles/Simulations", force = TRUE)
library(Simulations)
source("WrappersVblip1.R")

n = 5000
B = 96

dgps = lapply(1:B, FUN = function(x) get.dgp(n,4))

# detectCores()
# cl = makeCluster(2, type = "SOCK")
# registerDoSNOW(cl)
# clusterExport(cl,cl_export)

# debug(SL.stack1)
# debug(sim_cv)
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.hal")
# SL.library = list("SL.nnet", "glm.mainint", c("SL.hal", "screen.Main"))
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.glm.interaction", "SL.mean","SL.rpartPrune", "SL.earth", "SL.glmnet")
# SL.library = list("SL.nnet", "glm.mainint", "SL.mean", "SL.glm","SL.rpartPrune", "earthFull","SL.glmnet")
# SL.library = SL.libraryG = c("SL.glm", "SL.mean")



gform = formula("A~.")
Qform = formula("Y~A*(W1+W2+W3+W4)")
ALL = list()
for (a in 1:B) {
  S = sim_hal(data = dgps[[a]], gform = gform, Qform = Qform, 
              V = 10, single = FALSE, estimator = "single 1step", method = "method.NNLS", 
              gn = NULL, cvhal = TRUE, parallel = TRUE)
  ALL[[B]] = c(S, BV0 = dgp[[a]]$BV0, ATE0 = dgp[[a]]$ATE0)
}

save(ALL, dgps, file = "caseHalvsDelta5G.RData")

# res = lapply(ALL, FUN = function(x) x$res)
# results = data.matrix(data.frame(do.call(rbind, res)))
# 
# colnames(results)
# vapply(c(4, 10), FUN = function(x) {
#   mean(results[,x+1] <= results[,"ATE0"] & results[,x+2] >= results[,"ATE0"])
# }, FUN.VALUE = 1)
# 
# vapply(c(1, 7), FUN = function(x) {
#   mean(results[,x+1] <= results[,"BV0"] & 
#          results[,x+2] >= results[,"BV0"])
# }, FUN.VALUE = 1)
# 
# 
# # vapply(c(1, 41, 47), FUN = function(x) {
# #   results[,x+1] <= results[,"BV0"] & results[,x+2] >= results[,"BV0"]
# # }, FUN.VALUE = rep(1, nrow(results)))
# # grep("coef", colnames(results))
# #
# vapply(c(1, 7, 65, 77, 71), FUN = function(x) {
#   mean(results[,x] - results[, "BV0"])
# }, FUN.VALUE = 1)
# 
# vapply(c(4, 8, 68, 78, 74), FUN = function(x) {
#   mean(results[,x] - results[,"ATE0"])
# }, FUN.VALUE = 1)
# 
# BIAS = vapply(c(4, 8, 68, 78, 74), FUN = function(x) {
#   results[,x] - results[,"ATE0"]
# }, FUN.VALUE = rep(1, nrow(results)))
# 
# hist(BIAS[,1])
# 
# vapply(c(4, 68, 74), FUN = function(x) {
#   mean(results[,x+2] - results[,x+1])
# }, FUN.VALUE = 1)
# 
# vapply(c(1, 65, 71), FUN = function(x) {
#   mean(results[,x+2] - results[,x+1])
# }, FUN.VALUE = 1)
# #
# colMeans(results[,grep("coef", colnames(results))])