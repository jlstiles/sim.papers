
# function to check MSE
#' @export
perf=function(ests,truth){
  n=length(ests)
  var=((n-1)/n)*var(ests)
  bias=mean(ests)-truth
  mse=mean((ests-truth)^2)
  c(var=var,bias=bias,mse=mse)
}

#' @export
gendata=function(n,g0, Q0){
  W1 = runif(n,-3,3)
  W2=rnorm(n)
  W3=runif(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  data.frame(A,W1,W2,W3,W4,Y)
}

#' @export
Q0_trig1= function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 2*A * W1 + 20*cos(W1) * A - 3*W1 * sin(2*W2)+ cos(W1)
              -3*W2+4*A*(W2^2) +3*cos(W4)*A +A*W1^2- 2 * sin(W2)*W4 - 6*A* W3 * W4-3))
}

#' @export
Q0_trig =function (A, W1, W2, W3, W4)
{
  plogis(.14*(2* A  + 20*cos(W1) * A +cos(W1)-4*A*(W2^2) +3*cos(W4)*A +A*W1^2))
}

#' @export
Q0_1 = function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 2*A * W1 + 4*A*W3*W4+W2*W1+W3*W4+10*A*cos(W4)))
}


#' @export
g0_linear= function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4 - 0.15))
}

#' @export
g0_1 = function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.08 * W1^2*W2+.5*W1 + 0.49 * cos(W2)*W3 + 0.18 * W3^2 - 0.12 * sin(W4) - 0.15))
}

#' @export
sim_lr = function(n, g0, Q0, form) {
  
  simdata = gendata(n, g0, Q0)
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

#' @export
sim_hal = function(n, g0, Q0, HAL, SL.library, SL.libraryG) {
  
  simdata = gendata(n, g0, Q0)
  head(simdata)
  X=simdata
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newX = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  if (HAL==TRUE) {
    halresults <- hal(Y = simdata$Y,newX = newX,
                      X = X, family = binomial(),
                      verbose = FALSE, parallel = FALSE)
    
    Qk = halresults$pred[1:n]
    Q1k = halresults$pred[n+1:n]
    Q0k = halresults$pred[2*n+1:n]
    
    initest = var(Q1k-Q0k)
    
    gfit = glm(A~.,data = X, family = 'binomial')
    gk = predict(gfit, type = 'response')
    
    halresultsG <- hal(Y = simdata$A,newX = W,
                       X = W, family = binomial(),
                       verbose = FALSE, parallel = FALSE)
    gk_hal = halresultsG$pred[1:n]
  } else {
    data <- gendata(n)
    X = data
    X1 = X0 = X
    X0$A = 0
    X1$A = 1
    newdata = rbind(X,X1,X0)
    
    mainform = paste0(paste(colnames(data)[2:4],"+",collapse=""),colnames(data)[5])
    mainform
    squares = paste0(paste0("I(",colnames(data)[2:5]),"^2)")
    squares
    squares = paste0(paste(squares[1:3],"+",collapse=""),squares[4])
    squares
    mainsq = paste0(mainform,"+",squares)
    mainsq
    mainsq.int = paste0("Y~A*(",mainsq,")")
    mainsq.int = formula(mainsq.int) 
    mainsq.int
    
    newdata = model.matrix(mainsq.int,newdata)
    newdata = as.data.frame(newdata[,-1])
    colnames(newdata)[2:ncol(newdata)] = paste0("X",2:ncol(newdata))
    head(newdata)
    X = newdata[1:n,]
    X$Y = NULL
    
    # time = proc.time()
    Qfit=SuperLearner(data$Y,X,newX=newdata, family = binomial(),
                      SL.library=SL.library, 
                      id = NULL, verbose = FALSE, control = list(),
                      cvControl = list(V=10), obsWeights = NULL)
    
    # proc.time() - time
    W = X
    W$A = NULL
    
    gfit = SuperLearner(data$A,W,newW, family = binomial(),
                        SL.library=SL.libraryG, 
                        id = NULL, verbose = FALSE, control = list(),
                        cvControl = list(V=10), obsWeights = NULL)
    
    gk = gfit$library.predict[1:n,gfit$coef!=0] %*% gfit$coef[gfit$coef!=0]
    
    Qk = Qfit$library.predict[1:n,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    Q1k = Qfit$library.predict[(n+1):(2*n),Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    Q0k = Qfit$library.predict[(2*n+1):(3*n),Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    
  }
  
  
  initdata = data.frame(A = X$A, Y = simdata$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  initdataG = data.frame(A = X$A, Y = simdata$Y, gk = gk_hal, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  sigmait_infoG = gentmle2::gentmle(initdata=initdataG, params=list(param_sigmaATE), 
                                    submodel = submodel_logit, loss = loss_loglik,
                                    approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_infoG = gentmle2::gentmle(initdata=initdataG, params=list(param_sigmaATE), 
                                  submodel = submodel_logit, loss = loss_loglik,
                                  approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  results_glm <- glm(Y ~ A*(W1+W2+W3+W4),data = simdata,family=binomial())
  Qk = predict(results_glm, newdata = X, type = "response")
  Q1k = predict(results_glm,newdata=X1,type='response')
  Q0k = predict(results_glm,newdata=X0,type='response')
  initest_lr = var(Q1k - Q0k)
  initdata = data.frame(A = X$A, Y = simdata$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  
  sigmait_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                       submodel = submodel_logit, loss = loss_loglik,
                                       approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  steps = c(sigma_info$steps, sigma_infoG$steps,sigmait_info$steps, sigmait_infoG$steps,
            sigma_info_glm$steps,sigmait_info_glm$steps)
  converge = c(sigma_info$converge, sigma_infoG$converge,
               sigmait_info$converge, sigmait_infoG$converge,
               sigma_info_glm$converge, sigmait_info_glm$converge)
  
  ci_sig = ci_gentmle(sigma_info)[c(2,4,5)]
  ci_sigG = ci_gentmle(sigma_infoG)[c(2,4,5)]
  ci_sigit = ci_gentmle(sigmait_info)[c(2,4,5)]
  ci_sigitG = ci_gentmle(sigmait_infoG)[c(2,4,5)]
  ci_sig_glm = ci_gentmle(sigma_info_glm)[c(2,4,5)]
  ci_sigit_glm = ci_gentmle(sigma_info_glm)[c(2,4,5)]
  
  return(c(ci_sig,ci_sigG, ci_sigit, ci_sigitG, ci_sig_glm,
           ci_sigit_glm,initest, initest_lr,
           steps = steps, converge = converge))
}

