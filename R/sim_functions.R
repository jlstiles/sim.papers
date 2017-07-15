
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
sim_lr = function(n, g0, Q0, formQ, formG) {
  
  data = gendata(n, g0, Q0)
  
  X=data
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newX = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  results <- glm(formQ,data = data,family=binomial())
  Qk = predict(results, newdata = X, type = 'response')
  Q1k = predict(results, newdata = X1,type = 'response')
  Q0k = predict(results, newdata = X0,type = 'response')
  
  gfit = glm(formG,data = X, family = 'binomial')
  gk = predict(gfit, type = 'response')
  initest = var(Q1k - Q0k)
  initest_ATE = mean(Q1k - Q0k)
  
  initdata = data.frame(A = X$A, Y = data$Y, gk = gk, Qk = Qk, 
                        Q1k = Q1k, Q0k = Q0k)
  
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                       submodel = submodel_logit, loss = loss_loglik,
                                       approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE,param_sigmaATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                     simultaneous.inference = TRUE)
  
  ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2)
  
  steps = c(sigma_info$steps,sigmait_info$steps, simul_info$steps, ATE_info$steps)
  converge = c(sigma_info$converge, sigmait_info$converge,simul_info$converge, ATE_info$converge)
  
  ci_sig = ci_gentmle(sigma_info)[c(2,4,5)]
  ci_sigit = ci_gentmle(sigmait_info)[c(2,4,5)]
  ci_simul = ci_gentmle(simul_info)[2,c(2,4,5)]
  ci_simulATE = ci_gentmle(simul_info)[1,c(2,4,5)]
  ci_ATE = ci_gentmle(ATE_info)[c(2,4,5)]
  
  cis = c(ci_sig, ci_sigit, ci_simul, ci_simulATE, 
          ci_ATE)
  names(cis)[c(1,4,7,10,13)] =
    c("sig", "sigit","simul","simulATE","ATE")
  names(converge) = names(steps) = names(cis)[c(1:3,5)]
  results = c(cis, initest = initest, initest_ATE = initest_ATE,steps = steps, 
              converge = converge)
  
  return(results)
}

#' @export
sim_hal = function(n, g0, Q0, HAL, SL.library, SL.libraryG, method = "method.NNLS") {
  # SL.library = SL.libraryG = list("SL.glm","SL.mean")
  # method = "method.NNloglik"
  # HAL = FALSE
  # g0 = g0_1
  # Q0 = Q0_1
  # n=1000
  data = gendata(n, g0, Q0)
  # head(simdata)
  X=data
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newdata = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  results_glm <- glm(Y ~ A*(W1+W2+W3+W4),data = data,family=binomial())
  Qk_glm = predict(results_glm, newdata = X, type = "response")
  Q1k_glm = predict(results_glm,newdata=X1,type='response')
  Q0k_glm = predict(results_glm,newdata=X0,type='response')
  
  gfit_glm = glm(A~.,data = X, family = 'binomial')
  gk_glm = predict(gfit_glm, type = 'response')
  
  if (HAL==TRUE) {
    halresults <- hal(Y = data$Y,newX = newdata,
                      X = X, family = binomial(),
                      verbose = FALSE, parallel = FALSE)
    
    Qk = halresults$pred[1:n]
    Q1k = halresults$pred[n+1:n]
    Q0k = halresults$pred[2*n+1:n]
    
    initest = var(Q1k-Q0k)
    
    gfit = glm(A~.,data = X, family = 'binomial')
    gk = predict(gfit, type = 'response')
    
    halresultsG <- hal(Y = data$A,newX = W,
                       X = W, family = binomial(),
                       verbose = FALSE, parallel = FALSE)
    gk = halresultsG$pred[1:n]
    Qcoef = Gcoef = 0
  } else {
    
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
                      SL.library=SL.library, method=method,
                      id = NULL, verbose = FALSE, control = list(),
                      cvControl = list(V=10), obsWeights = NULL)
    
    
    gfit = SuperLearner(data$A,W,newX = W, family = binomial(),
                        SL.library=SL.libraryG,method = method, 
                        id = NULL, verbose = FALSE, control = list(),
                        cvControl = list(V=10), obsWeights = NULL)
    
    if (length(gfit$coef[gfit$coef!=0])==1){
      gk = gfit$library.predict[1:n,gfit$coef!=0]
      } else {
      gk = gfit$library.predict[1:n,gfit$coef!=0] %*% gfit$coef[gfit$coef!=0]
      }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Qk = Qfit$library.predict[1:n,Qfit$coef!=0]
    } else {
      Qk = Qfit$library.predict[1:n,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Q1k = Qfit$library.predict[n+1:n,Qfit$coef!=0]
    } else {
      Q1k = Qfit$library.predict[n+1:n,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Q0k = Qfit$library.predict[2*n+1:n,Qfit$coef!=0]
    } else {
      Q0k = Qfit$library.predict[2*n+1:n,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    Qcoef = Qfit$coef
    Gcoef = gfit$coef
    
    Qrisk = Qfit$cvRisk
    grisk = gfit$cvRisk
    
  }
  
  initest = var(Q1k - Q0k)
  initest_ATE = mean(Q1k - Q0k)
  initdata = data.frame(A = X$A, Y = data$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2)
  
  simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                 simultaneous.inference = TRUE)
  
  simuljl_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "line", max_iter = 100,g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
  
  simuljer_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
  
  ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                    submodel = submodel_logit, loss = loss_loglik,
                                    approach = "full", max_iter = 100,g.trunc = 1e-2)
  
  initest_lr = var(Q1k_glm - Q0k_glm)
  initest_lr_ATE = mean(Q1k_glm - Q0k_glm)
  
  initdata = data.frame(A = X$A, Y = data$Y, gk = gk_glm, Qk = Qk_glm, 
                        Q1k = Q1k_glm, Q0k = Q0k_glm)
  
  sigmait_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                       submodel = submodel_logit, loss = loss_loglik,
                                       approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  simul_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_ATE,param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 10000, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
  
  ATE_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "full", max_iter = 100, g.trunc = 1e-2)
  
  steps = c(sigma_info$steps, sigmait_info$steps, simul_info$steps, simuljl_info$steps,
            simuljer_info$steps, sigma_info_glm$steps,sigmait_info_glm$steps, 
            simul_info_glm$steps,simul_info$steps, ATE_info$steps, simul_info_glm$steps
            ,ATE_info_glm$steps)
  converge = c(sigma_info$converge, sigmait_info$converge,simul_info$converge, 
               simuljl_info$converge, simuljer_info$converge, sigma_info_glm$converge,
               sigmait_info_glm$converge,simul_info_glm$converge,simul_info$converge, 
               ATE_info$converge, simul_info_glm$converge, ATE_info_glm$converge)

  ci_sig = ci_gentmle(sigma_info)[c(2,4,5)]
  ci_sigit = ci_gentmle(sigmait_info)[c(2,4,5)]
  ci_simul = ci_gentmle(simul_info)[2,c(2,4,5)]
  ci_simuljl = ci_gentmle(simuljl_info)[2,c(2,4,5)]
  ci_simuljer = ci_gentmle(simuljer_info)[2,c(2,4,5)]
  ci_sig_glm = ci_gentmle(sigma_info_glm)[c(2,4,5)]
  ci_sigit_glm = ci_gentmle(sigmait_info_glm)[c(2,4,5)]
  ci_simul_glm = ci_gentmle(simul_info_glm)[2,c(2,4,5)]
  
  ci_simulATE = ci_gentmle(simul_info)[1,c(2,4,5)]
  ci_ATE = ci_gentmle(ATE_info)[c(2,4,5)]  
  ci_simulATE_glm = ci_gentmle(simul_info_glm)[1,c(2,4,5)]
  ci_ATE_glm = ci_gentmle(ATE_info_glm)[c(2,4,5)]
  
  cis = c(ci_sig,ci_sigit, ci_simul, ci_simuljl, ci_simuljer,
    ci_sig_glm, ci_sigit_glm, ci_simul_glm, ci_simulATE,ci_ATE, ci_simulATE_glm, 
    ci_ATE_glm)
  names(converge) = names(steps) = names(cis)[c(1,4,7,10,13,16,19,22,25,28,31,34)]=
    c("sig", "sigit", "simul", "simul_line", "simul_full","sig_glm", 
      "sigit_glm","simul_glm","simulATE","ATE","simulATE_glm","ATE_glm")
  results = c(cis, initest = initest, initest_lr = initest_lr,initest_ATE = initest_ATE, 
              initest_lr_ATE = initest_lr_ATE, steps = steps, converge = converge, Qcoef, 
              Gcoef, Qrisk, grisk)
  # results
  return(results)
}

#' @export
sim_single = function(n, g0, Q0, form) {
  
  X = gendata(n, g0, Q0)
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newX = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  results_glm <- glm(form,data = simdata,family = binomial())
  Qk = predict(results_glm, newdata = X, type = 'response')
  Q1k = predict(results_glm, newdata = X1,type = 'response')
  Q0k = predict(results_glm, newdata = X0,type = 'response')
  initest = var(Q1k - Q0k)
  
  gfit = glm(A~.,data = X, family = binomial())
  gk = predict(gfit, type = 'response')
  
  initdata = data.frame(A = X$A, Y = simdata$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  ci_tmle = ci_gentmle(sigma_info)[c(2,4,5)]
  return(c(ci_tmle, initest))
}

#' @export
SL.stack1 = function(Y, X, A, W, newdata, method, SL.library, SL.libraryG, ...) {
  # 
  # X = X
  # Y = data$Y
  # A = data$A
  # W = X
  # W$A = NULL
  # newdata = newdata
  # method = "method.NNloglik"
  folds = make_folds(n, V=10)
  stack = lapply(folds, FUN = function(x) {
    # x=folds[[5]]
    tr = x$training_set
    val = x$validation_set
    n=length(tr)
    nv = length(val)
    
    Y = Y[tr]
    X = X[tr,]
    newtr = c(val, (n+val),(2*n+val))
    newdata = newdata[newtr,]
    Qfit=SuperLearner(Y,X,newX=newdata, family = binomial(),
                      SL.library=SL.library, method=method,
                      id = NULL, verbose = FALSE, control = list(),
                      cvControl = list(V=10), obsWeights = NULL)
    
    A = A[tr]
    W1 = W[tr,]
    newW = W[val,]
    gfit = SuperLearner(Y=A,X=W1,newX = newW, family = binomial(),
                        SL.library=SL.libraryG,method = method, 
                        id = NULL, verbose = FALSE, control = list(),
                        cvControl = list(V=10), obsWeights = NULL)
    
    
    if (length(gfit$coef[gfit$coef!=0])==1){
      gk = gfit$library.predict[1:nv,gfit$coef!=0]
    } else {
      gk = gfit$library.predict[1:nv,gfit$coef!=0] %*% gfit$coef[gfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Qk = Qfit$library.predict[1:nv,Qfit$coef!=0]
    } else {
      Qk = Qfit$library.predict[1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Q1k = Qfit$library.predict[nv+1:nv,Qfit$coef!=0]
    } else {
      Q1k = Qfit$library.predict[nv+1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Q0k = Qfit$library.predict[2*nv+1:nv,Qfit$coef!=0]
    } else {
      Q0k = Qfit$library.predict[2*nv+1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    Qcoef = Qfit$coef
    Gcoef = gfit$coef
    
    Qrisk = Qfit$cvRisk
    Grisk = gfit$cvRisk
    
    return(list(Qk = Qk, Q0k = Q0k, Q1k = Q1k, gk = gk, Qcoef = Qcoef, Gcoef = Gcoef,
                Qrisk = Qrisk, Grisk = Grisk, inds = x$validation_set))
  })
  
  if (length(SL.library)==1) {
    Qcoef = mean(Qcoef_mat)
    Qrisk = mean(Qrisk_mat)
  } else {
    Qcoef = rowMeans(Qcoef_mat)
    Qrisk = rowMeans(Qrisk_mat)
  }
  
  if (length(SL.libraryG)==1) {
    Gcoef = mean(Gcoef_mat)
    Grisk = mean(Grisk_mat)
  } else {
    Gcoef = rowMeans(Gcoef_mat)
    Grisk = rowMeans(Grisk_mat)
  }
  
  inds = unlist(lapply(stack, FUN = function(x) x$inds))
  Y = Y[inds]
  A = A[inds]
  
  initdata = data.frame(Y=Y,A=A,Qk=Qk,Q1k=Q1k,Q0k=Q0k,gk=gk)
  return(list(initdata = initdata, Qcoef = Qcoef, Gcoef = Gcoef, Qrisk = Qrisk,
              Grisk = Grisk, inds = inds))
}  
# 
#' @export
SL.stack = function(Y, X, A, W, newdata, method, SL.library, SL.libraryG, mc.cores = 1, ...) {
  # 
  # X = X
  # Y = data$Y
  # A = data$A
  # W = X
  # W$A = NULL
  # newdata = newdata
  # method = "method.NNloglik"
  folds = make_folds(n, V=10)
  stack = mclapply(folds, FUN = function(x) {
    # x=folds[[5]]
    tr = x$training_set
    val = x$validation_set
    n=length(tr)
    nv = length(val)
    
    Y = Y[tr]
    X = X[tr,]
    newtr = c(val, (n+val),(2*n+val))
    newdata = newdata[newtr,]
    Qfit=SuperLearner(Y,X,newX=newdata, family = binomial(),
                      SL.library=SL.library, method=method,
                      id = NULL, verbose = FALSE, control = list(),
                      cvControl = list(V=10), obsWeights = NULL)
    
    A = A[tr]
    W1 = W[tr,]
    newW = W[val,]
    gfit = SuperLearner(Y=A,X=W1,newX = newW, family = binomial(),
                        SL.library=SL.libraryG,method = method, 
                        id = NULL, verbose = FALSE, control = list(),
                        cvControl = list(V=10), obsWeights = NULL)
    
    
    if (length(gfit$coef[gfit$coef!=0])==1){
      gk = gfit$library.predict[1:nv,gfit$coef!=0]
    } else {
      gk = gfit$library.predict[1:nv,gfit$coef!=0] %*% gfit$coef[gfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Qk = Qfit$library.predict[1:nv,Qfit$coef!=0]
    } else {
      Qk = Qfit$library.predict[1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Q1k = Qfit$library.predict[nv+1:nv,Qfit$coef!=0]
    } else {
      Q1k = Qfit$library.predict[nv+1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1){
      Q0k = Qfit$library.predict[2*nv+1:nv,Qfit$coef!=0]
    } else {
      Q0k = Qfit$library.predict[2*nv+1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    Qcoef = Qfit$coef
    Gcoef = gfit$coef
    
    Qrisk = Qfit$cvRisk
    Grisk = gfit$cvRisk
    
    return(list(Qk = Qk, Q0k = Q0k, Q1k = Q1k, gk = gk, Qcoef = Qcoef, Gcoef = Gcoef,
                Qrisk = Qrisk, Grisk = Grisk, inds = x$validation_set))
  }, mc.cores = mc.cores)
  
  Qk = unlist(lapply(stack, FUN = function(x) x$Qk))
  Q1k = unlist(lapply(stack, FUN = function(x) x$Q1k))
  Q0k = unlist(lapply(stack, FUN = function(x) x$Q0k))
  gk = unlist(lapply(stack, FUN = function(x) x$gk))
  
  Qcoef_mat = vapply(stack, FUN = function(x) x$Qcoef, 
                     FUN.VALUE = rep(1,length(stack[[1]]$Qcoef)))
  Qrisk_mat = vapply(stack, FUN = function(x) x$Qrisk,
         FUN.VALUE = rep(1,length(stack[[1]]$Qrisk)))
  
  Gcoef_mat = vapply(stack, FUN = function(x) x$Gcoef, 
                     FUN.VALUE = rep(1,length(stack[[1]]$Gcoef)))
  Grisk_mat = vapply(stack, FUN = function(x) x$Grisk,
                      FUN.VALUE = rep(1,length(stack[[1]]$Gcoef)))
  
  if (length(SL.library)==1) {
    Qcoef = mean(Qcoef_mat)
    Qrisk = mean(Qrisk_mat)
  } else {
    Qcoef = rowMeans(Qcoef_mat)
    Qrisk = rowMeans(Qrisk_mat)
  }
  
  if (length(SL.libraryG)==1) {
    Gcoef = mean(Gcoef_mat)
    Grisk = mean(Grisk_mat)
  } else {
    Gcoef = rowMeans(Gcoef_mat)
    Grisk = rowMeans(Grisk_mat)
  }

  inds = unlist(lapply(stack, FUN = function(x) x$inds))
  Y = Y[inds]
  A = A[inds]
  
  initdata = data.frame(Y=Y,A=A,Qk=Qk,Q1k=Q1k,Q0k=Q0k,gk=gk)
  return(list(initdata = initdata, Qcoef = Qcoef, Gcoef = Gcoef, Qrisk = Qrisk,
              Grisk = Grisk, inds = inds))
}  

#' @export
sim_cv = function(n, g0, Q0, SL.library, SL.libraryG, method = "method.NNLS") {
  
  # n=1000
  # g0 = g0_linear
  # Q0 = Q0_trig
  # SL.library = SL.libraryG = c("SL.mean", "SL.glm")
  # method = "method.NNloglik"
  
  data = gendata(n, g0, Q0)
  
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
  Y = data$Y
  A = data$A
  W = X[,2:5]
  
  # time = proc.time()
  stack = SL.stack1(Y=Y, X=X, A=A, W=W, newdata=newdata, method=method, 
                      SL.library=SL.library, SL.libraryG=SL.libraryG)
  # proc.time() - time

  initdata = stack$initdata   
  initest = with(initdata,var(Q1k - Q0k))
  initest_ATE = with(initdata, mean(Q1k - Q0k))
  
  
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2)
  
  simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                 simultaneous.inference = TRUE)
  
  simuljl_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "line", max_iter = 100,g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
  
  simuljer_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                    submodel = submodel_logit, loss = loss_loglik,
                                    approach = "full", max_iter = 100,g.trunc = 1e-2,
                                    simultaneous.inference = TRUE)
  
  ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                               submodel = submodel_logit, loss = loss_loglik,
                               approach = "full", max_iter = 100,g.trunc = 1e-2)
  
  results_glm <- glm(Y ~ A*(X2+X3+X4+X5),data = X[,1:5],family=binomial())
  Q = predict(results_glm, newdata = newdata[,1:5], type = "response")
  Qk_glm = Q[1:n]
  Q1k_glm = Q[n+1:n]
  Q0k_glm = Q[2*n+1:n]
  gfit_glm = glm(A~.,data = X[,1:5], family = 'binomial')
  gk_glm = predict(gfit_glm, type = 'response')
  initest_lr = var(Q1k_glm - Q0k_glm)
  initest_lr_ATE = mean(Q1k_glm - Q0k_glm)
  
  initdata = data.frame(A = X$A, Y = data$Y, gk = gk_glm, Qk = Qk_glm, 
                        Q1k = Q1k_glm, Q0k = Q0k_glm)
  
  sigmait_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                       submodel = submodel_logit, loss = loss_loglik,
                                       approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  simul_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_ATE,param_sigmaATE), 
                                     submodel = submodel_logit, loss = loss_loglik,
                                     approach = "full", max_iter = 10000, g.trunc = 1e-2,
                                     simultaneous.inference = TRUE)
  
  ATE_info_glm = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2)
  
  steps = c(sigma_info$steps, sigmait_info$steps, simul_info$steps, simuljl_info$steps,
            simuljer_info$steps, sigma_info_glm$steps,sigmait_info_glm$steps, 
            simul_info_glm$steps,simul_info$steps, ATE_info$steps, simul_info_glm$steps
            ,ATE_info_glm$steps)
  converge = c(sigma_info$converge, sigmait_info$converge,simul_info$converge, 
               simuljl_info$converge, simuljer_info$converge, sigma_info_glm$converge,
               sigmait_info_glm$converge,simul_info_glm$converge,simul_info$converge, 
               ATE_info$converge, simul_info_glm$converge, ATE_info_glm$converge)
  
  ci_sig = ci_gentmle(sigma_info)[c(2,4,5)]
  ci_sigit = ci_gentmle(sigmait_info)[c(2,4,5)]
  ci_simul = ci_gentmle(simul_info)[2,c(2,4,5)]
  ci_simuljl = ci_gentmle(simuljl_info)[2,c(2,4,5)]
  ci_simuljer = ci_gentmle(simuljer_info)[2,c(2,4,5)]
  ci_sig_glm = ci_gentmle(sigma_info_glm)[c(2,4,5)]
  ci_sigit_glm = ci_gentmle(sigmait_info_glm)[c(2,4,5)]
  ci_simul_glm = ci_gentmle(simul_info_glm)[2,c(2,4,5)]
  
  ci_simulATE = ci_gentmle(simul_info)[1,c(2,4,5)]
  ci_ATE = ci_gentmle(ATE_info)[c(2,4,5)]  
  ci_simulATE_glm = ci_gentmle(simul_info_glm)[1,c(2,4,5)]
  ci_ATE_glm = ci_gentmle(ATE_info_glm)[c(2,4,5)]
  
  cis = c(ci_sig,ci_sigit, ci_simul, ci_simuljl, ci_simuljer,
          ci_sig_glm, ci_sigit_glm, ci_simul_glm, ci_simulATE,ci_ATE, ci_simulATE_glm, 
          ci_ATE_glm)
  names(converge) = names(steps) = names(cis)[c(1,4,7,10,13,16,19,22,25,28,31,34)]=
    c("sig", "sigit", "simul", "simul_line", "simul_full","sig_glm", 
      "sigit_glm","simul_glm","simulATE","ATE","simulATE_glm","ATE_glm")
  results = c(cis, initest = initest, initest_lr = initest_lr,initest_ATE = initest_ATE, 
              initest_lr_ATE = initest_lr_ATE, steps = steps, converge = converge, 
            Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
            Grisk = stack$Grisk)
  # results
  return(results)
}

# pp = sim_hal(100, g0 = g0_1, Q0 = Q0_trig, HAL = TRUE, SL.library=NULL, SL.libraryG=NULL)
# pp[c(1,4,7,10,13,16,19,22,23)]
