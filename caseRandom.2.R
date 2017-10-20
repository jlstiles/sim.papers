case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jeremyrcoyle/gentmle2")
# devtools::install_github("jlstiles/Simulations", force = TRUE)
library(Simulations)
source("WrappersVblip1.R")

SL.library = SL.libraryD2
SL.libraryG = SL.libraryGD2

n = 1000
B = 100

dgps = lapply(1:B, FUN = function(x) {
  get.dgp(n = n, d = 2, pos = 0.01, minATE = -2, minBV = 0.03, depth = 2, maxterms = 2, minterms = 1, 
          mininters = 1, num.binaries = 0, force.confounding = TRUE)
})

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)

# debug(SL.stack1)
# debug(sim_cv)
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.hal")
# SL.library = list("SL.nnet", "glm.mainint", c("SL.hal", "screen.Main"))
# SL.libraryG = c("SL.glm", "SL.nnet", "SL.glm.interaction", "SL.mean","SL.rpartPrune", "SL.earth", "SL.glmnet")
# SL.library = list("SL.nnet", "glm.mainint", "SL.mean", "SL.glm","SL.rpartPrune", "earthFull","SL.glmnet")
# SL.library = c("glm.mainint", "nnetMain")
# SL.libraryG = c("SL.glm", "SL.nnet")

gform = formula("A~.")
Qform = formula("Y~A*(W1+W2)")
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = g0_linear, Q0 = Q0_trig1, SL.library = SL.library,
                    SL.libraryG = SL.libraryG, method = "method.NNloglik", cv = TRUE, V = 10, SL = 10L,
                    gform = gform, Qform = Qform, estimator = c("single 1step"), dgp = dgps[[i]]
            )}


save(ALL, dgps, file = "caseRandomD2.2.RData")

