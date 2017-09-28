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
ALL=foreach(i=1:10,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = NULL, Q0 = Q0_trig1, SL.library = SL.library,
                    SL.libraryG = SL.libraryG, method = "method.NNLS", cv = TRUE, V = 10, SL = 10L,
                    gform = gform, Qform = Qform, estimator = c("single 1step"), dgp = NULL, gn = NULL
            )}

save(ALL, dgps, file = "caseRandom.1.RData")