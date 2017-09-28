case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jeremyrcoyle/gentmle2")
# devtools::install_github("jlstiles/Simulations", force = TRUE)
library(Simulations)
source("WrappersVblip1.R")

SL.library = SL.library1
SL.libraryG = SL.libraryG

n = 1000
B = 100

dgps = lapply(1:B, FUN = function(x) get.dgp(1000,4))

length(dgps)

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)

# debug(SL.stack1)
# debug(sim_cv)
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.hal")
# SL.library = list("SL.nnet", "glm.mainint", c("SL.hal", "screen.Main"))
# SL.libraryG = c("SL.glm")
# SL.library = list("SL.nnet", "glm.mainint", "SL.mean")

gform = formula("A~.")
Qform = formula("Y~A*(W1+W2+W3+W4)")
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = NULL, Q0 = Q0_trig1, SL.library = SL.library,
                    SL.libraryG = SL.libraryG, method = "method.NNLS", cv = TRUE, V = 10, SL = 10L,
                    gform = gform, Qform = Qform, estimator = c("single 1step"), dgp = dgps[[i]], gn = NULL
            )}

save(ALL, dgps, file = "caseRandom.RData")

# results = data.matrix(data.frame(do.call(rbind, ALL)))
# 
# res = lapply(ALL, FUN = function(x) x$res)
# results = data.matrix(data.frame(do.call(rbind, res)))
# 
# vapply(c(4, 24, 30), FUN = function(x) {
#   mean(results[,x+1] <= results[,"ATE0"] & results[,x+2] >= results[,"ATE0"])
# }, FUN.VALUE = 1)
# 
# vapply(c(1, 21, 27), FUN = function(x) {
#   mean(results[,x+1] <= results[,"BV0"] & results[,x+2] >= results[,"BV0"])
# }, FUN.VALUE = 1)
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
