# basic example with very simple SuperLearner library
SL.library = c("SL.glm", "SL.mean")
SL.libraryG = c("SL.glm", "SL.mean")
n=1000
result = sim_cv(n, g0 = g0_1, Q0 = Q0_1, SL.library=SL.library,
       SL.libraryG=SL.libraryG,
       method = "method.NNloglik",cv=TRUE)
          
names(result)
result
