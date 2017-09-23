case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jlstiles/Simulations")
library(Simulations)
source("WrappersVblip1.R")

SL.library = SL.library3
SL.libraryG = "SL.glm"

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)
n = 1000
B = 100

# SL.library = c("SL.glm", "SL.mean")
# SL.libraryG = "SL.glm"
g0 = g0_linear
Q0 = Q0_trig1
# debug(SL.stack1)
# debug(sim_cv)
gform = formula("A~.")
Qform = formula("Y~A*(W1+W2+W3+W4)")
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = g0, Q0 = Q0, SL.library = SL.library,
                    SL.libraryG = SL.libraryG, method = "method.NNloglik", cv = TRUE, V = 10, SL = 10L, 
                    gform = gform, Qform = Qform, estimator = c("single 1step")
            )}
results = data.matrix(data.frame(do.call(rbind, ALL)))
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
save(ALL, file = "case2aSL3.RData")


