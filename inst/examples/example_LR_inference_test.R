
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4))

# get a randomly drawn dataframe under the specified model
simmie = function(n) {
  data = gendata(n, g0_linear, Q0_linear)
  # well-specified model
  Qform = formula("Y ~ W1 + W2 + A*(W3 + W4)")
  
  # specifying the covariates, treatment and outcome
  W = data[,2:5]
  A = data$A
  Y = data$Y
  
  # should cover each truth 95 percent of the time.
  info = LR.TSM(W=W,A=A,Y=Y,Qform=Qform, setA = 1, alpha = .05)
  truth = get.truth(g0_linear, Q0_linear)
  cover = (info[5] < truth[2]) & (info[6] > truth[2])
  return(cover)
}

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl,cl_export)

B = 100
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {simmie(1000)}
mean(unlist(ALL))

# should cover each truth 95 percent of the time and both truths
# simultaneously 95 percent of the time for the simultaneous CI's
info1 = LR.inference(W=W,A=A,Y=Y,Qform=Qform, alpha = .05, 
                    simultaneous.inference = TRUE)
info1
