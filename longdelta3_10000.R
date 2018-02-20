library(Simulations)
library(simcausal)
options(simcausal.verbose=FALSE)
t.end = 7


D <- DAG.empty()
D <- D +
  node("L2", t = 0, distr = "rbern",
       prob = 0.45) +
  # node("L1", t = 0, distr = "rbern",
  #      prob = ifelse(L2[0] == 1, 0.5, 0.1)) +
  node("A1", t = 0, distr = "rbern",
       prob = plogis(-1 + .8*L2[0]))

# define variables based on past values in time
D <- D +
  node("Y", t = 1:t.end, distr = "rbern",
       prob =
         plogis(-2.5 +  4 * L2[t-1] + 0.05 * sum(I(L2[0:(t-1)] == rep(0, t)))),
       EFU = TRUE) +
  node("L2", t = 1:t.end, distr = "rbern",
       prob =
         ifelse(A1[t-1] == 1, 0.1,
                ifelse(L2[t-1] == 1, 0.9, min(1, 0.1 + t / 16)))) +
  node("A1", t = 1:t.end, distr = "rbern",
       prob = ifelse(A1[t-1] == 1, .4,
                     ifelse(L2[t] == 0, 0.3, 0.6)))
lDAG <- set.DAG(D)
Ddyn <- lDAG

act_theta0 <-node("A1", t = 0:2, distr = "rbern",
                  prob = (t==0)*1 + (t==1)*0 + (t==2)*1)

Ddyn <- Ddyn + action("A1_th0", nodes = act_theta0)

Ynodes = as.character(paste0("Y_",1:7))
Anodes = as.character(paste0("A1_",0:6))
Lnodes = as.character(paste0("L2_",0:6))
# specify the formulas
formula0 = formula("Y_1 ~ A1_0*L2_0")
formula1 = formula("Y_2 ~ L2_0*L2_1*A1_0*A1_1")
formula2 = formula("Y_3 ~ L2_2*A1_2*L2_1*A1_0*A1_1*L2_0")
formula3 = formula("Y_4 ~ L2_3*A1_3")
formula4 = formula("Y_5 ~ L2_4*A1_4")
formula5 = formula("Y_6 ~ L2_5*A1_5")
formula6 = formula("Y_7 ~ L2_6*A1_6")
formulas = list(formula0, formula1, formula2, formula3, formula4,
                formula5, formula6)

Qform <- c(Y_1="Q.kplus1 ~ A1_0*L2_0", Y_2="Q.kplus1 ~ L2_1*A1_0*A1_1*L2_0",
           Y_3="Q.kplus1 ~ L2_2*A1_2*L2_1*A1_0*A1_1*L2_0")
gform <- c(A1_0="A1_0 ~ L2_0",A1_1="A1_1 ~ L2_1*A1_0",
           A1_2="A1_2 ~ A1_1*L2_2")

sim.longTSM = function(n, dag, gform, Qform, formulas, setA, T_end, 
                       Lnodes, Anodes, Ynodes)
{
  
  OdatL = sim(Ddyn, n = n)
  OdatL$ID = NULL
  data = OdatL
  data_ltmle = data
  for (t in 1:(T_end-1)) {
    data_ltmle[data_ltmle[,Ynodes[t]]==1,Ynodes[t+1]] = 1 
  }
  
  nombre = Ynodes[T_end]
  Yend = grep(nombre, colnames(data))
  
  res = ltmle(data=data_ltmle[,1:Yend], Anodes=Anodes[1:T_end], Lnodes = Lnodes[1:T_end], 
              Ynodes=Ynodes[1:T_end],survivalOutcome = TRUE, abar = setA[1:T_end], 
              Qform = Qform[1:T_end], gform = gform[1:T_end], 
              gbounds = c(0.000001,1),deterministic.g.function = NULL,  
              estimate.time = TRUE, gcomp = TRUE, iptw.only = FALSE, stratify = FALSE,
              deterministic.Q.function = NULL,variance.method = "ic", 
              observation.weights = NULL, id = NULL)
  
  TSMinfo = long.TSM(data = data, Ynodes = Ynodes[1:T_end], Anodes = Anodes[1:T_end], 
                     formulas = formulas[1:T_end], setA = setA[1:T_end])
  
  TSMinfo$CI
  sd(TSMinfo$IC)*sqrt(n-1)/n
  summary(res)[[1]]$std.dev
  sd(res$IC$iptw)*sqrt(n-1)/n
  
  c(summary(res)[[1]]$estimate, summary(res)[[1]]$CI)
  
  CIs = c(c(summary(res)[[1]]$estimate, summary(res)[[1]]$CI),summary(res)[[1]]$std.dev,
          TSMinfo$CI, sd(TSMinfo$IC)*sqrt(n-1)/n,sd(res$IC$iptw)*sqrt(n-1)/n)
  names(CIs)[c(2:3,6:7)] = c("left", "right")
  names(CIs)[c(1,4,5,8,9)] = c("gcomp", "SE gcomp","LRdelta","SE LR", "SE iptw")
  return(CIs)
}

cl = makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

# run this on a 24 core node
B=1000
n=10000
setA = c(1,0,1,1,1,1,0)
T_end = 3

ALL=foreach(i=1:B,.packages=c("Simulations","simcausal","parallel"),
            .errorhandling = "remove")%dopar%
            {sim.longTSM(n=n, dag=Ddyn, gform=gform, Qform=Qform, 
                         formulas =formulas, setA = setA, T_end = 3,
                         Lnodes = Lnodes, Anodes = Anodes, Ynodes = Ynodes)}

ALL3_10000=ALL

Ddyn <- Ddyn + action("A1_th0", nodes = act_theta0)

Ddyn1 <- set.targetE(Ddyn, outcome = "Y", t = 3, param = "A1_th0")
psi0 = eval.target(Ddyn1, n = 1000000)$res
psi0

res = do.call(rbind,ALL3_10000)
coverage = cov.check(res,psi0,c(1,5))
coverage
save(formulas, Qform, gform, Ddyn, coverage, psi0,setA,
     Ynodes, Anodes, Lnodes,  ALL3_10000, T_end, file = "longdelta3_10000.RData")
