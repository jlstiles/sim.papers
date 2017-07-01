
# function to check MSE
perf=function(ests,truth){
  n=length(ests)
  var=((n-1)/n)*var(ests)
  bias=mean(ests)-truth
  mse=mean((ests-truth)^2)
  c(var=var,bias=bias,mse=mse)
}

gendata=function(n,g0, Q0){
  W1 = runif(n,-3,3)
  W2=rnorm(n)
  W3=runif(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  data.frame(A,W1,W2,W3,W4,Y)
}


Q0_trig1= function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 2*A * W1 + 20*cos(W1) * A - 3*W1 * sin(2*W2)+ cos(W1)
              -3*W2+4*A*(W2^2) +3*cos(W4)*A +A*W1^2- 2 * sin(W2)*W4 - 6*A* W3 * W4-3))
}

Q0_trig =function (A, W1, W2, W3, W4)
{
  plogis(.14*(2* A  + 20*cos(W1) * A +cos(W1)-4*A*(W2^2) +3*cos(W4)*A +A*W1^2))
}

Q0_1 = function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 2*A * W1 + 4*A*W3*W4+W2*W1+W3*W4+10*A*cos(W4)))
}



g0_linear= function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4 - 0.15))
}

g0_1 = function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.08 * W1^2*W2+.5*W1 + 0.49 * cos(W2)*W3 + 0.18 * W3^2 - 0.12 * sin(W4) - 0.15))
}

sim_lr = function(n, g0, Q0) {
  
  simdata = gendata(n, g0, Q0, form)
  head(simdata)
  X=simdata
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newX = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  results_glm <- glm(form,data = simdata,family=binomial())
  Qk = predict(results_glm,newdata = X, type = 'response')
  Q1k = predict(results_glm,newdata=X1,type='response')
  Q0k = predict(results_glm,newdata=X0,type='response')
  initest = var(Q1k-Q0k)
  return(initest)
}
