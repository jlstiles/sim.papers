
# devtools::install_github("jlstiles/Simulations")
library(Simulations)
source("wrappers_ex.R")

g0 = function (W1, W2) {
  plogis(.4*(-0.4 * W1*W2 + 0.63 * W2^2 -.66*cos(W1) - 0.25))
}

Q0 = function (A, W1, W2) {
  plogis(0.2 * W1*W2 + 0.1 * W2^2 - .8*A*(cos(W1) + .5*A*W1*W2^2) - 0.35)
}

gendata.fcn = function (n, g0, Q0) 
{
  W1 = runif(n, -3, 3)
  W2 = rnorm(n)
  A = rbinom(n, 1, g0(W1, W2))
  Y = rbinom(n, 1, Q0(A, W1, W2))
  data.frame(A, W1, W2, Y)
}

# pop = gendata.fcn(1e6, g0, Q0)
# pscores = with(pop, g0(W1, W2))
# hist(pscores, 200)
# 
# blips = with(pop, Q0(1, W1, W2) - Q0(0, W1, W2))
# hist(blips, 200)
# var(blips)
# Q1s = with(pop, Q0(1, W1, W2))
# Q0s = with(pop, Q0(0, W1, W2))
# Qs = with(pop, Q0(A, W1, W2))
# 
# hist(Q1s, 200)
# hist(Q0s, 200)
# hist(Qs, 200)


SL.libraryD2 = list("nnetMain","nnetMain1","glm.mainint", 
                    "earth_2d","SL.glm.interaction", "xgboost_2d","SL.mean","SL.hal")
SL.libraryGD2 = list("nnetMain","nnetMain1", "earth_2d","SL.glm.interaction", "xgboost_2dG","SL.mean","SL.hal")

cl_export = c("nnetMain","nnetMain1","glm.mainint", "earth_2d","xgboost_2d","xgboost_2dG", "SL.hal")

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)
n = 1000
B = 100


# debug(SL.stack1)
# debug(sim_cv)
gform = formula("A~.")
Qform = formula("Y~A*(W1+W2)")
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"), 
            .errorhandling = "remove") %dopar%
            {sim_cv(n, g0 = g0, Q0 = Q0, SL.library = SL.libraryD2,
                    SL.libraryG = SL.libraryD2G, method = "method.NNloglik", cv = TRUE, V = 10, SL = 10L, 
                    gform = gform, Qform = Qform, estimator = c("single 1step"), gendata.fcn = gendata.fcn)
            }
# results = data.matrix(data.frame(do.call(rbind, ALL)))

save(ALL, file = "case_d24.RData")


