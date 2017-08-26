# Simulations
Read the [Blip Variance Paper, full version](https://github.com/jlstiles/Simulations/blob/master/blip-variance-article.pdf) 

To install gentmle and the other tools necessary to obtain results in in the
Blip Variance paper:

devtools::install_github("jlstiles/Simulations")

The following file has instructions for how to generate all of the results in the full version of the Blip Variance paper.  I will note to the reader in the code comments, how to proceed with parallelizing so the process does not take so much time as well as issues with memory one should consider when running ensemble methods as applied here.

[SET UP: run this first](#setup)

























<a name="setup"></a>
**Set Up**
Set the source file paths as below for source_paper, WrappersVblip1.R, Wrappers_ex.R (for the example), and the wcgs.dta file.  Always run this chunk first.

```R
case = "setup"
resultsGotten = TRUE
source_file = "~/Dropbox/Jonathan/Simulations/source_paper.R"
source(source_file)

devtools::install_github("jlstiles/Simulations")
library(Simulations)
source("~/Dropbox/Jonathan/Simulations/WrappersVblip1.R")
# source('/home/jlstiles/R/WrappersVblip1.R')
```

```R
# case "lr_case2a"

case = "LRcase2a"
resultsGotten = FALSE
B=3
n=1000
no.cores = detectCores()
source(source_file)

# head of compiled simulation results, all is computed from this data.frame 
head(results_LRcase2a) 

# figure in the paper
ggoverLRcase2a

```

```R
case = "LRcase2b"
resultsGotten = FALSE
B=3
n=1000
no.cores = detectCores()
source(source_file)

# compiled simulation results, all is computed from this data.frame 
head(results_LRcase2b)

# figure in the paper
gg_LRcase2b
```

```R
case = "HALcase2a"
resultsGotten = FALSE
B=3
n=100
no.cores = detectCores()
source(source_file)

results_HALcase2a
gg_HALcase2a
stargazer(coverage_HALcase2a, summary = FALSE)
stargazer(performance.sig_HALcase2a, summary = FALSE)
```

```R
# THIS CHUNK FOR HAL tmle vs LR tmle, case 2b
case = "HALcase2b"
resultsGotten = FALSE
B=3
n=100
no.cores = detectCores()
source(source_file)

results_HALcase2b = results
gg_HALcase2b
stargazer(coverage_HALcase2b, summary = FALSE)
stargazer(performance.sig_HALcase2b, summary = FALSE)
```


```R
# case2a
case = "case2a"
resultsGotten = TRUE
B=n=1000
no.cores = detectCores()
source(source_file)
 
gg_ATEcase2a
gg_BVcase2a
stargazer(performance.ate_case2a,summary=FALSE,digits=6)
stargazer(performance.sig_case2a,summary=FALSE,digits=6)
stargazer(coverage_case2a,summary=FALSE,digits=3)
stargazer(SL_results_case2a,summary=FALSE,digits=5)
```

```R
# THIS CHUNK IS FOR CASE 2b using SuperLearner Library 1, SL1
case = "case2bSL1"
resultsGotten = FALSE
B=3
n=1000
no.cores = detectCores()
source(source_file)

results_case2bSL1 = results
gg_ATEcase2bSL1
gg_BVcase2bSL1
stargazer(performance.ate_case2bSL1,summary=FALSE,digits=6)
stargazer(performance.sig_case2bSL1,summary=FALSE,digits=6)
stargazer(coverage_case2bSL1,summary=FALSE,digits=3)
stargazer(SL_results_case2bSL1,summary=FALSE,digits=5)
```

```R
#####
##### other sim--CV-tMLE includes the asshole learners

# RUN ten times with B  = 100 on 24 cores each
g0 = g0_linear
Q0 = Q0_trig
testdata=gendata(1000000, g0=g0, Q0 = Q0)
blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
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
results = results_case2bCVSL2
case = "case2bCVSL2"
source(source_file)
gg_ATEcase2bCVSL2
gg_BVcase2bCVSL2
stargazer(performance.ate_case2bCVSL2,summary=FALSE,digits=6)
stargazer(performance.sig_case2bCVSL2,summary=FALSE,digits=6)
stargazer(coverage_case2bCVSL2,summary=FALSE,digits=3)
stargazer(SL_results_case2bCVSL2,summary=FALSE,digits=5)
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
# THIS CHUNK IS FOR CASE 2b using SuperLearner Library 2, SL2,
# which overfit and caused skewing necessitating CV-TMLE
case = "case2bSL2"
resultsGotten = FALSE
B=3
n=1000
no.cores = detectCores()
source(source_file)

# save for later compilation
results_of = results

ggover_ATEcase2bSL2
ggover_BVcase2bSL2
stargazer(performance.ate_case2bSL2,summary=FALSE,digits=6)
stargazer(performance.sig_case2bSL2,summary=FALSE,digits=6)
stargazer(coverage_case2bSL2,summary=FALSE,digits=3)
stargazer(SL_results_case2bSL2,summary=FALSE,digits=5)

results[order(results_case2bSL2[,"Qcoef.rangerFull_screen.Main"]),1]

df = data.frame(x=results_case2bSL2[,"Qcoef.xgbFull_All"], 
                y=results_case2bSL2[,1])
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
#####
##### case 4 
case = "case4"
resultsGotten = FALSE
B=3
n=100
no.cores = detectCores()
source(source_file)

results_case4 = results
# ggover_ATEcase4
ggover_BVcase4
stargazer(MSE_cov_case4,summary=FALSE,digits=6)
stargazer(SL_results_case4,summary=FALSE,digits=5)
```

```R
#####
##### case 3
case = "case3"
resultsGotten = FALSE
B=3
n=1000
no.cores = detectCores()
source(source_file)

ggover_ATEcase3
ggover_BVcase3
stargazer(performance.ate_case3,summary=FALSE,digits=6)
stargazer(performance.sig_case3,summary=FALSE,digits=6)
stargazer(coverage_case3,summary=FALSE,digits=3)
stargazer(SL_results_case3,summary=FALSE,digits=5)
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

save(gg_pscoresWell, ml250, ml500, ml1000, gg_wells,
     file = "~/Dropbox/Jonathan/Simulations/results/wells.RData")
```

```R
# noise simulations 
# set case = "noise" or case = "noise_neg" for negative bias
case = "noise"
# case = "noise_neg"
no.cores  = 4
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

save(gg_pscoresWell, ml, ml1, gg_bias, gg_mse, gg_coverage,L,
     file = "~/Dropbox/Jonathan/Simulations/results/noise.RData")

```

```R
# noise simulations 
# set case = "noise" or case = "noise_neg" for negative bias
# case = "noise"
case = "noise_neg"
no.cores  = 4
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

save(gg_pscoresWell, ml, ml1, gg_bias, gg_mse, gg_coverage,L,
     file = "~/Dropbox/Jonathan/Simulations/results/noiseneg.RData")
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
