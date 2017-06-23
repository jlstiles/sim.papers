library(knitr)
library(dplyr)
library(reshape)
library(reshape2)
library(plyr)
library(ggplot2)
library(Hmisc)
library(stargazer)
library(xtable)
library(ltmle)
library(SuperLearner)
library(gridExtra)
library(cowplot)
library(tmle)
library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(nnet)
library(caret)
library(e1071)
library(gam)
library(rpart)
library(randomForest)
library(gentmle2)

g0= function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4 - 0.15))
}

gendata=function(n){
  W1 = runif(n,-3,3)
  # W1= rnorm(n)
  # W1=rnorm(n)
  W2=rbinom(n,1,.5)
  W3=rnorm(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  data.frame(A,W1,W2,W3,W4,Y)
}

Q0 =function (A, W1, W2, W3, W4)
  {
    plogis(.2*(.1*A+2*A*W1-10*A*W2+3*A*W3^2+W1+W2+.4*W3+.3*W4))
  }
  
# making sure true barQ is not idiotic
truth = gendata(1000000)
Q_true = with(truth, Q0(A,W1,W2,W3,W4))
hist(Q_true)

Q1_true = with(truth, Q0(1,W1,W2,W3,W4))
Q0_true = with(truth, Q0(0,W1,W2,W3,W4))
blip_true = Q1_true-Q0_true
hist(blip_true)
var0 = var(blip_true)
ATE0 = mean(blip_true)
ATE0
var0

# Correlating normals so the variance is not increased when subtracting. This
# makes compatible noise on both barQ_1 and barQ_0 and barQ
alpha = .5
z1 = rnorm(1000000)
z2 = alpha*z1+(1-alpha^2)^.5*rnorm(1000000)
cor(z1,z2)
sd(z2)
sd(z1)
var(z1-z2)
a = rbinom(1000000,1,.5)
Z = a*z1 + (1-a)*z2
var(Z)

# add noise to Q_true
truth = gendata(1000000)
noise_analysis = function(n, rate, truth) {
  # rate = .33
  # # V=1
  # # truth = gendata(10000)
  # n
  Q1_true = with(truth, Q0(1,W1,W2,W3,W4))
  Q0_true = with(truth, Q0(0,W1,W2,W3,W4))
  Q_true = with(truth,Q0(A,W1,W2,W3,W4))
  # g_true = with(truth, g0(W1,W2,W3,W4))
  
  # biasG = function(W1,W2,W3,W4,n,rate)  n^-rate*5*(.5+.2*W1+1*W2*.8*W3-2*W4)
  # hist(with(truth,biasG(W1,W2,W3,W4,n=n)))
  biasQ = function(A,W1,W2,W3,W4,n,rate)  -n^-rate*.8*(+.2+1.5*A+.2*W1+1*W2+5*A*W3^2+ 1*W4)
  # hist(with(truth,biasQ(A,W1,W2,W3,W4,n=n,rate=rate)))
  sdQ = function(A,W1,W2,W3,W4,n,rate) (n^-rate)*.8*(abs(3.5+.5*W1+.15*W2+.33*W3*W4-W4))
  # hist(with(truth,sdQ(A,W1,W2,W3,W4,n=n,rate=rate)))
  # sdG = function(W1,W2,W3,W4,V=0,n,rate) (n^-rate)*1*(abs(1.5-.1*W2+.5*W4+.25*W3-.23*W4*W1-W4)/300+V)
  # hist(with(truth,sdG(W1,W2,W3,W4,V=V,n=n)))
  
  # sdQ = function(A,W1,W2,W3,W4,V=0,n,rate) 0
  # sdG = function(W1,W2,W3,W4,V=0,n,rate) 0
  noise.1 = function(data,n,rate) rnorm(nrow(data), with(data, biasQ(1,W1,W2,W3,W4,n=n,rate=rate)),
                                 with(data, sdQ(1,W1,W2,W3,W4,n=n,rate=rate)))
  # noise on barQ(0,W) is correlated per draw so not any more variant
  noise.0 = function(data,noise,n,rate) .5*noise + 
    sqrt(.75)*rnorm(nrow(data), with(data, biasQ(0,W1,W2,W3,W4,n=n,rate=rate)),
                    with(data, sdQ(0,W1,W2,W3,W4,n=n,rate=rate)))
  
  noise = function(data, noise_1, noise_0) with(data, A*noise_1+(1-A)*noise_0)
  # noiseG = function(data,V,n,rate) rnorm(nrow(data), with(data, biasG(W1,W2,W3,W4,n=n,rate=rate)),
                                                 # with(data, sdG(W1,W2,W3,W4,V=V,n=n,rate=rate)))
  noise_1 = noise.1(truth,n,rate)
  noise_0 = noise.0(truth,noise_1,n,rate)
  noise_A = noise(truth, noise_1, noise_0)
  # noise_G = noiseG(truth,V,n,rate)
  var(noise_1-noise_0)
  var(noise_A)
  Q1_test = plogis(qlogis(Q1_true) + noise_1)
  Q0_test = plogis(qlogis(Q0_true) + noise_0)
  Q_test = plogis(qlogis(Q_true) + noise_A)
  # g_test = plogis(qlogis(g_true) + noise_G)
  
  blip_test = Q1_test - Q0_test
  blip_true = Q1_true - Q0_true
  
  L2_blip = sqrt(mean((blip_test - blip_true)^2))
  L2_Q = sqrt(mean((Q_true - Q_test)^2))
  # L2_G = sqrt(mean((g_true - g_test)^2))
  
  hist(Q_test)
  hist(Q_true)
  hist(blip_test)
  hist(blip_true)
  # hist(g_test)
  ate_bias = mean(blip_test) - ATE0
  var_bias = var(blip_test) - var0
  
  title = paste("sample size ",n, "variance bias = ",round(var_bias,5))
  df = data.frame(blip = c(blip_test,blip_true), type = c(rep("test",1e6),rep("true",1e6)))
  gg_testvstrue = ggplot(df,aes(x=blip,color=type))+geom_density()++ggtitle(title)+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size=8,angle=315))  
  # var_bias
  # ate_bias
  
  results = list(L2_blip = n^.25*L2_blip, L2_Q = n^.25*L2_Q,
              ate_bias = ate_bias, var_bias = var_bias,plot = gg_testvstrue)
  return(results)
}

# undebug(noise_analysis)
# seq(100,10000,1000)
noise_look = vapply(rep(250,10),FUN = function(x) noise_analysis(n=x,rate=1/3,truth=truth), 
                    FUN.VALUE = c(1,1,1,1))

noise_look
rowMeans(noise_look)

# truth = gendata(1000000)
vapply(1:10,FUN = function(x) noise_analysis(10000,1,1/3,truth), FUN.VALUE = rep(1,5))

plot(seq(250,100000,10000),noise_look[1,])
plot(seq(250,100000,10000),noise_look[2,])
plot(seq(250,100000,10000),noise_look[3,])
plot(seq(250,100000,10000),noise_look[4,])


# We now set up a simulation with a similar architecture to create bias

# get truths
truth = gendata(1000000)
Q1_true = with(truth, Q0(1,W1,W2,W3,W4))
Q0_true = with(truth, Q0(0,W1,W2,W3,W4))
Q_true = with(truth,Q0(A,W1,W2,W3,W4))
g_true = with(truth, g0(W1,W2,W3,W4))

# biasG = function(W1,W2,W3,W4,n,rate)  n^-rate*5*(.5+.2*W1+1*W2*.8*W3-2*W4)
# hist(with(truth,biasG(W1,W2,W3,W4,n=n)))
biasQ = function(A,W1,W2,W3,W4,n,rate)  -n^-rate*.8*(+.2+1.5*A+.2*W1+1*W2+5*A*W3^2+ 1*W4)
# hist(with(truth,biasQ(A,W1,W2,W3,W4,n=n,rate=rate)))
sdQ = function(A,W1,W2,W3,W4,n,rate) (n^-rate)*.8*(abs(3.5+.5*W1+.15*W2+.33*W3*W4-W4))
# hist(with(truth,sdQ(A,W1,W2,W3,W4,n=n,rate=rate)))
# sdG = function(W1,W2,W3,W4,V=0,n,rate) (n^-rate)*3*(abs(1.5-.1*W2+.5*W4+.25*W3-.23*W4*W1-W4)/300+V)
# hist(with(truth,sdG(W1,W2,W3,W4,V=V,n=n)))

noise.1 = function(data,n,rate) rnorm(nrow(data), with(data, biasQ(1,W1,W2,W3,W4,n=n,rate=rate)),
                                        with(data, sdQ(1,W1,W2,W3,W4,n=n,rate=rate)))
# noise on barQ(0,W) is correlated per draw so not any more variant
noise.0 = function(data,noise,n,rate) .5*noise + 
  sqrt(.75)*rnorm(nrow(data), with(data, biasQ(0,W1,W2,W3,W4,n=n,rate=rate)),
                  with(data, sdQ(0,W1,W2,W3,W4,n=n,rate=rate)))

noise = function(data, noise_1, noise_0) with(data, A*noise_1+(1-A)*noise_0)
# noiseG = function(data,V,n,rate) rnorm(nrow(data), with(data, biasG(W1,W2,W3,W4,n=n,rate=rate)),
                                       # with(data, sdG(W1,W2,W3,W4,V=V,n=n,rate=rate)))
# Simulation function
simBlipvar = function(n,rate){
  # tack on the noise and use as an initial estimate
  data = gendata(n)
  noise_1 = noise.1(data,n,rate)
  noise_0 = noise.0(data, noise_1,n,rate)
  noise_A = noise(data, noise_1, noise_0)
  # noise_G = noiseG(data,V,n,rate)
  
  Q1k = plogis(qlogis(with(data,Q0(1,W1,W2,W3,W4)))+noise_1)
  Q0k = plogis(qlogis(with(data,Q0(0,W1,W2,W3,W4)))+noise_0)
  Qk = plogis(qlogis(with(data,Q0(A,W1,W2,W3,W4)))+noise_A)
  # gk = plogis(qlogis(with(data,g0(W1,W2,W3,W4)))+noise_G)
  gk = with(data,g0(W1,W2,W3,W4))
  
  initdata = data.frame(Qk=Qk,Q1k=Q1k,Q0k=Q0k,gk=gk,A=data$A,Y=data$Y)
  
  ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                               submodel = submodel_logit, loss = loss_loglik,
                               approach = "line", max_iter = 100,g.trunc = 1e-2)
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  steps = c(ATE_steps = ATE_info$steps, sigmait_steps = sigmait_info$steps,
            sigma_info = sigma_info$steps)
  
  ATE_ci = gentmle2::ci_gentmle(ATE_info)[c(2,4,5)]
  sigmait_ci = gentmle2::ci_gentmle(sigmait_info)[c(2,4,5)]
  sigma_ci = gentmle2::ci_gentmle(sigma_info)[c(2,4,5)]
  
  converges = c(sigmait = sigmait_info$converge,sigma = sigma_info$converge,
                ATE = ATE_info$converge)
  
  initest = var(Q1k-Q0k)
  return(c(sigmait_ci = sigmait_ci, ci_sigma = sigma_ci, ATE_ci = ATE_ci,
           sigma_init = initest, ATE_init = ATE_info$initests,
           steps = steps, converges = converges))
}

# check the function works with for each
detectCores()
cl = makeCluster(4, type = "SOCK")
registerDoSNOW(cl)


B=100
n=10000

rate = 1/3
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm"))%dopar%
{simBlipvar(n=n,rate=rate)}

getRes(ALL,100)


