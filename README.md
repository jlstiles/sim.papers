# Simulations
Read the [Blip Variance Paper, full version](https://github.com/jlstiles/Simulations/blob/master/blip-variance-article.pdf) 

To install gentmle and the other tools necessary to obtain results in in the
Blip Variance paper:

devtools::install_github("jlstiles/Simulations")

The following knitr file has instructions for how to generate the results, including plots and charts, in the paper.  I will note to the reader in the code comments, how to proceed with parallelizing so the process does not take so much time as well as issues with memory one should consider when running ensemble methods as applied here.

[SET UP: run this first](#setup)

[section 6.1.2 manufactured noise simulations](#section6.1.2)

[section 6.1.3 manufactured noise simulations](#section6.1.3)

[section 6.3 well specified precaution on skewing](#section6.3)

[section 6.4.1 narrow model failures](#section6.4.1)

[section 6.4.2 narrow model failures](#section6.4.2)

well-spec pscore, mispec. outcome model

[section 6.4.3 case 2a](#section6.4.3)

[section 6.4.4 case2b SuperLearner lib 1](#section6.4.42bSL1)

[section 6.4.4 case2b SuperLearner lib 2 with CV-TMLE](#section6.4.42bCVSL2)

[section 6.4.4 case2b SuperLearner lib 2 with TMLE, overfit hurts without CV](#section6.4.4 2bSL2)

mispecified treatment and outcome models

[section 6.4.4 combining cases as in paper](#section6.4.4 2bcombo)

[section 6.5.1](#section6.5.1)

[section 6.5.2](#section6.5.2)

[section 7 real data example](#section7)

<a name="setup"></a>
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

<a name="section6.4.1"></a>
```R
# section 6.4.1
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

<a name="section6.4.2"></a>
```R
# section 6.4.2
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

<a name="section6.4.3"></a>
```R
# section 6.4.3
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

<a name="section6.4.42bSL1"></a>
```R
# section 6.4.4 using SuperLearner Library 1
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

<a name="section6.4.42bCVSL2"></a>
```R
# section 6.4.4 using SuperLearner Library 2 but with CV-TMLE

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

<a name="section6.4.42bSL2"></a> 
```R
# section 6.4.4 using SuperLearner Library 2 without CV-tmle so 
# it has a skewed sampling dist

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

<a name="section6.4.42bcombo"></a>   
```R
# section 6.4.4
# Have the following variables named as follows 
# results_case2bSL1 = results for case 2b with SL.library1
# results_case2bSL2 = results for case 2b with SL.library2, the overfitting lib.
# results_case2bCVSL2 = results for case 2b with SL.library2, but using CV-TMLE

case = "combo_case2b"
source(source_file)

ggover_cvadvert = ggover
stargazer(MSE_cov)
```       

<a name="section6.5.1"></a>  
```R
# section 6.5.1
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

<a name="section6.5.2"></a>  
```R
# section 6.5.2
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

<a name="section6.3"></a>  
```R
##### section 6.3
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

<a name="section6.1.2"></a>  
```R
# noise simulations section 6.1.2
# set case = "noise" or case = "noise_neg" for negative bias
case = "noise"
# case = "noise_neg"
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

save(gg_pscoresWell, ml, ml1, gg_bias, gg_mse, gg_coverage,L,
     file = "~/Dropbox/Jonathan/Simulations/results/noise.RData")

```

<a name="section6.1.3"></a>  
```R
# noise simulations section 6.1.3
# set case = "noise" or case = "noise_neg" for negative bias
# case = "noise"
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

save(gg_pscoresWell, ml, ml1, gg_bias, gg_mse, gg_coverage,L,
     file = "~/Dropbox/Jonathan/Simulations/results/noiseneg.RData")
```

```R
# To obtain CIs via delta method for logistic regression plug-in estimates
# step 1:choose case
LRcases = list("CI_LRcase2a", "CI_LRcase2b", "CI_LRcase3", "CI_LRcase4")
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

<a name="section7"></a>  
```R
# section 7
case = "example"
library(foreign)
library(Simulations)

source('/home/jlstiles/R/wrappers_ex.R')
wcgs = read.dta("/home/jlstiles/R/wcgs.dta")

# This will set the SL libraries.  We transform the design matrix
# to include squares, main terms and interactions and screen from there
source(source_file)

# This will take a couple of days as below. We limited the cores to 2 
# due to RAM issues but youc an specify the number of cores in the function below
stack = SL.stack(Y = Y, X = X, A = A, W = W, newdata = newdata, 
                 method = "method.NNloglik",
                 SL.library = SL.library, 
                 SL.libraryG = SL.libraryG, V=10, mc.cores = 2)

# perform the targeting step.  See package examples for more info
info = gentmle(stack$initdata, params = list(param_ATE,param_sigmaATE),
               approach = "full", max_iter = 100, simultaneous.inference = TRUE)
info$tmleests
info$initests
results = rbind(info$initests, info$tmleests)

# getting estimate for standard dev of blip using delta method
psi_sd = info$tmleests[2]^(.5)
IC_sd = .5*psi_sd^(-1)*info$Dstar[,2]

# log scaling blip variance due to left bound of CI in neg range
psi_log = log(info$tmleests[2])
IC_log = 1/info$tmleests[2]*info$Dstar[,2]

# correlation matrix the same whether using log scale, stand dev of blip
# or blip variance since all of these ICs are constants times each other
IC_ate = info$Dstar[,1]
IC = data.frame(IC_sd = IC_sd, IC_ate = IC_ate)
Sig = cor(IC)

# getting simultaneous zscore
library(mvtnorm)
z = rmvnorm(1000000, c(0,0), Sig)
n = length(IC[,2])
zscore = quantile(apply(z,1,FUN = function(x) max(abs(x))),.95)
zscore

# This ci for blip variance automatically generated by the gentmle
ci = ci_gentmle(info)

# other cis
ci_simul_sd = c(psi_sd, psi_sd-zscore*sd(IC_sd)*sqrt(n-1)/n, 
                psi_sd+zscore*sd(IC_sd)*sqrt(n-1)/n)

ci_simul_log = exp(c(psi_log,psi_log - zscore*sd(IC_log)*sqrt(n-1)/n,
                   psi_log + zscore*sd(IC_log)*sqrt(n-1)/n))

# compiling cis
ci_ate_sd_log = rbind(ci[1,c(2,4,5)], ci_simul_sd, ci_simul_log, ci[2,c(2,4,5)])
rownames(ci_ate_sd_log) = c("ate", "blip st. dev","blip var log-scaled","blip var")

library(stargazer)
stargazer(ci_ate_sd_log, summary= FALSE, digits =5)

# steps to convergence
info$steps
info$converge

# superlearner coefficients 
stack$Qcoef
stack$Gcoef


```