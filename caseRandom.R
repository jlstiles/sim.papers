case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jlstiles/Simulations")
library(Simulations)
source("WrappersVblip1.R")

SL.library = SL.library1
SL.libraryG = SL.libraryG

n = 1000
B = 100

dgps = lapply(1:B, FUN = function(x) {
  info = get.info(1000,4,TRUE)
})

info = lapply(dgps, FUN = function(x) {
  list(DF = x$DF, BV0 = x$BV0, ATE0 = x$ATE0)
})

# detectCores()
# cl = makeCluster(detectCores(), type = "SOCK")
# registerDoSNOW(cl)
# clusterExport(cl,cl_export)
# library(doParallel)
# cl <- makePSOCKcluster(24)
# registerDoParallel(cl)
# clusterExport(cl,cl_export)

registerDoParallel(24)
# debug(SL.stack1)
# debug(sim_cv)
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.hal")
# SL.library = list("SL.nnet", "glm.mainint", c("SL.hal", "screen.Main"))
# SL.libraryG = c("SL.glm", "SL.nnet")
# SL.library = list("SL.nnet", "glm.mainint")

gform = formula("A~.")
Qform = formula("Y~A*(W1+W2+W3+W4)")
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = NULL, Q0 = Q0, SL.library = SL.library,
                    SL.libraryG = SL.libraryG, method = "method.NNloglik", cv = TRUE, V = 10, SL = 10L, 
                    gform = gform, Qform = Qform, estimator = c("single 1step"), dgp = info[[i]]
            )}

# results = data.matrix(data.frame(do.call(rbind, ALL)))
# results
# ALL = list()
# for (it in 1:6) {
#   ALL[[it]] = foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
#                      .errorhandling = "remove")%dopar%
#                      {sim_cv(n, g0 = g0, Q0 = Q0, SL.library = SL.library, 
#                              SL.libraryG = SL.libraryG[c(1:3,5:7)[1:it]], method = "method.NNLS", cv = TRUE, V = 2, SL = 2L, single = TRUE
#                      )}}
# 
# lapply(ALL, length)
save(ALL, dgps, file = "caseRandom.RData")

