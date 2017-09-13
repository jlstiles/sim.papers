
case = "setup"
source_file = "source_paper.R"
source(source_file)

# devtools::install_github("jlstiles/Simulations")
library(Simulations)
source("WrappersVblip1.R")

SL.library = SL.library1
SL.libraryG = SL.libraryG

cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)
n=1000
B=100
g0 = g0_1
Q0 = Q0_2
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                    SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=TRUE
            )}
results = data.matrix(data.frame(do.call(rbind, ALL)))

save(results, ALL, file = "case3.3.RData")