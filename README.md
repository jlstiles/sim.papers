# Simulations
Simulation tools

devtools::install_github("jlstiles/Simulations")

The following knitr file has instructions for how to generate the results, including plots and charts, in the paper.  I will note to the reader in the code comments, how to proceed with parallelizing so the process does not take so much time as well as issues with memory one should consider when running ensemble methods as applied here.

```R 
case = "setup"
resultsGotten = TRUE
source_file = "~/Dropbox/Jonathan/Simulations/source_paper.R"
source(source_file)

devtools::install_github("jlstiles/Simulations")
library(Simulations)
source("~/Dropbox/Jonathan/Simulations/CopyOfR/WrappersVblip1.R")
# source('/home/jlstiles/R/WrappersVblip1.R')
```

```R
# case "lr_case2a"

case = "lr_case2a"
resultsGotten = "TRUE"
B=n=1000
no.cores = detectCores()
source(source_file)

results
ggoverLR1
```

```R
case = "lr_case2b"
resultsGotten = "TRUE"
B=n=1000
no.cores = detectCores()
source(source_file)
library(stargazer)
results
ggoverLR2
```


<<eval = FALSE, echo = FALSE>>==
case = "hal_case2a"
resultsGotten = "TRUE"
B=n=1000
no.cores = detectCores()
source(source_file)

results
ggoverHAL2a
stargazer(coverage, summary = FALSE)
stargazer(performance.sig, summary = FALSE)
```

<<eval = FALSE, echo = FALSE>>==
# THIS CHUNK FOR HAL tmle vs LR tmle, case 2b
case = "hal_case2b"
resultsGotten = "TRUE"
B=n=1000
no.cores = detectCores()
source(source_file)

results
ggoverHAL2b
stargazer(coverage, summary = FALSE)
stargazer(performance.sig, summary = FALSE)
```


```R
# case2a
case = "case2a"
resultsGotten = TRUE
B=n=1000
no.cores = detectCores()
source(source_file)
 
ggover_ATEcase2a
ggover_BVcase2a
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)
```

```R
# THIS CHUNK IS FOR CASE 2b using SuperLearner Library 1, SL1
case = "case2b"
resultsGotten = TRUE
B=n=1000
no.cores = detectCores()
source(source_file)

results_SL1 = results
ggover_ATEcase2b
ggover_BVcase2b
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)
```

```R
# THIS CHUNK IS FOR CASE 2b using SuperLearner Library 2, SL2,
# which overfit and caused skewing necessitating CV-TMLE
case = "case2b_OF"
resultsGotten = TRUE
B=n=1000
no.cores = detectCores()
source(source_file)

# save for later compilation
results_of = results

ggover_ATEcase2bOF
ggover_BVcase2bOF
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)

results[order(results[,"Qcoef.rangerFull_screen.Main"]),1]

df = data.frame(x=results[,"Qcoef.xgbFull_All"], y=results[,1])
assholeForest = ggplot(df, aes(x=x,y=y))+geom_point()+
  labs(x="Random Forest Superlearner Coeff",y="blip variance estimate")+
  geom_vline(xintercept = var0, color = "green")
caps = paste0("Truth is green line\n",
              "The bigger the random Forest coefficient the worse outliers due to\n",
              "overfitting the outcome regression")

assholeForest=ggdraw(add_sub(assholeForest,caps, 
                       x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                       vpadding = grid::unit(1, "lines"), 
                       fontfamily = "", fontface = "plain",
                       colour = "black", size = 10, angle = 0, lineheight = 0.9))
assholeForest
```

```R
# has asshole learners but sample size was 2000 so the assholes didn't win any coeff,
# hence performance is awesome--included in previous chunck but here has more detail

# run this script for B=500 two times with 12 cores each 
# if RAM is 64G
g0 = g0_linear
Q0 = Q0_trig
testdata=gendata(1000000, g0=g0, Q0 = Q0)
blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
propensity = with(testdata, g0(W1,W2,W3,W4))
ATE0 = mean(blip_true)
var0 = var(blip_true)

SL.library = SL.library1
SL.libraryG = list("SL.glm")

cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)

ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                    SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
            )}

results = data.matrix(data.frame(do.call(rbind, ALL)))
results1 = data.matrix(data.frame(do.call(rbind, ALL)))
results = rbind(results, results1)
# save for later compilation
results_2G = results

#####
####
####

# rbind all the results files into one fill called results then run
case = "case2b_2G"
source(source_file) 

ggover_ATEcase2b2G
ggover_BVcase2b2G
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)

results[order(results[,"Qcoef.rangerFull_screen.Main"]),1]

hist(results[results[,"Qcoef.rangerFull_screen.Main"]==0,1])
```


```R
#####
##### other sim--CV-tMLE includes the asshole learners

# RUN ten times with B  = 100 on 24 cores each
g0 = g0_linear
Q0 = Q0_trig
testdata=gendata(1000000, g0=g0, Q0 = Q0)
blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
propensity = with(testdata, g0(W1,W2,W3,W4))
ATE0 = mean(blip_true)
var0 = var(blip_true)

SL.library = SL.library2
SL.libraryG = list("SL.glm")

cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)

B=1
n=1000
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                    SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
            )}

results = data.matrix(data.frame(do.call(rbind, ALL)))

# Once you have compiled 10 of these results via rbind, into one results 
# data.frame do the following for later combining of results

###
###
###
###

# rbind all 10 results data.frames into one data.frame named results
results = results_cv
case = "case2b_CV"
source(source_file)
ggover_ATEcase2bCV
ggover_BVcase2bCV
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)
  ```
 
```R
# For this chunk, I show the reader how to combine the results for case 2b
# Have the following variables named as follows 
# results = results for case 2b with SL.library1
# results_of = results for case 2b with SL.library2, the overfitting lib.
# results_cv = results for case 2b with SL.library2, but using CV-TMLE
# results_2G = results for case 2b with SL.library2, but sample size n = 2000
# 
case = "combo_case2b"
source(source_file)

ggover_cvadvert = ggover
stargazer(MSE_cov)
```       
 
  
```R
#####
##### case 4 
case = "case4"
resultsGotten = "TRUE"
B=n=1000
no.cores = detectCores()
source(source_file)

results
ggover_ATEcase4
ggover_BVcase4
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)
```

```R
#####
##### case 3
case = "case3"
resultsGotten = TRUE
B=n=1000
no.cores = detectCores()
source(source_file)

ggover_ATEcase3
ggover_BVcase3
stargazer(performance.ate,summary=FALSE,digits=6)
stargazer(performance.sig,summary=FALSE,digits=6)
stargazer(coverage,summary=FALSE,digits=3)
stargazer(SL_results,summary=FALSE,digits=5)
```

```R
#####
##### wells
case = "wells"
source(source_file)

# to print the three arranged plots, each with four plots
nums = c(1,2,3,4)

# sample size 250
ml250=marrangeGrob(gg_wells[c(nums)],ncol=2,nrow=2, widths = c(3.5,3.5),
                heights=rep(c(1,1)))

# sample size 500
ml500=marrangeGrob(gg_wells[c(nums+16)],ncol=2,nrow=2, widths = c(3.5,3.5),
                heights=rep(c(1,1)))

# sample size 1000
ml1000=marrangeGrob(gg_wells[c(nums+32)],ncol=2,nrow=2, widths = c(3.5,3.5),
                heights=rep(c(1,1)))

# propensity score density
gg_pscoresWell
# see sim_well_small.R and writeup_wells.Rnw
```

```R
# noise simulations 
# set case = "noise" or case = "noise_neg" for negative bias
case = "noise"
case = "noise_neg"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)

plotdf  = data.frame(n = sizes,coverage = res_noise[,5])
gg_coverage = ggplot(plotdf, aes(x = n, y = coverage)) + geom_point() +
  geom_hline(yintercept = .95, color = "green")

gg_bias = ggplot(plotdf_biasmse, aes(x = n, y = bias, color = type)) + 
  geom_line() + geom_hline(yintercept = 0, color = "green")

gg_mse = ggplot(plotdf_biasmse, aes(x = n, y = mse, color = type)) + 
  geom_line() + geom_hline(yintercept = 0, color = "green")

# plots of sampling dists
ml

# plots of blip dists true overlayed with blip dists tmle
ml1

```

```R
# noise simulations 
# set case = "noise" or case = "noise_neg" for negative bias
case = "noise"
case = "noise_neg"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)

plotdf  = data.frame(n = sizes,coverage = res_noise[,5])
gg_coverage = ggplot(plotdf, aes(x = n, y = coverage)) + geom_point() +
  geom_hline(yintercept = .95, color = "green")

gg_bias = ggplot(plotdf_biasmse, aes(x = n, y = bias, color = type)) + 
  geom_line() + geom_hline(yintercept = 0, color = "green")

gg_mse = ggplot(plotdf_biasmse, aes(x = n, y = mse, color = type)) + 
  geom_line() + geom_hline(yintercept = 0, color = "green")

# plots of sampling dists
ml

# plots of blip dists true overlayed with blip dists tmle
ml1

```

```R
# To obtain CIs via delta method for logistic regression plug-in estimates
# step 1:choose case
LRcases = list("LRcase2a", "LRcase2b", "LRcase3", "LRcase4")
case = "LRcase2a"
# number of simulations
B=1000
# sample size of each simulation
n=1000
source(source_file)
coverageLR
# you may rename it according to case
assign(paste0("coverage",case), coverageLR)
```