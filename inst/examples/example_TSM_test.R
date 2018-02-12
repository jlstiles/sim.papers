
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4)+W3+W4)
big = gendata(1e6,g0_linear, Q0_linear)
setA = 1
truth1 = mean(with(big, Q0_linear(setA, W1,W2,W3,W4)))
truth1
truth0 = mean(with(big, Q0_linear(0, W1,W2,W3,W4)))
truth0
# get a randomly drawn dataframe under the specified model
simmie = function(n, truth) {
  # n=1000
  data = gendata(n, g0_linear, Q0_linear)
  # well-specified model
  Qform = formula("Y ~ W1 + W2 + A*(W3 + W4)")
  
  # specifying the covariates, treatment and outcome
  W = data[,2:5]
  A = data$A
  Y = data$Y
  
  # should cover each truth 95 percent of the time.
  info = LR.TSM(W=W,A=A,Y=Y,Qform=Qform, setA = 0, alpha = .05)
  
  cover = (info$CI[2] < truth) & (info$CI[3] > truth)
  # cover
  # info$CI
  
  return(cover)
}

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)

B = 100
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {simmie(1000, truth0)}
results =  mean(unlist(ALL))
results

# should cover each truth 95 percent of the time and both truths
# simultaneously 95 percent of the time for the simultaneous CI's

