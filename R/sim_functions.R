library(boot)

#' @title SL.stack
#' @description Function that runs fully nested SuperLearner 
#' cross validated estimates 
#' on V-folds.  Only supports binary treatment.
#' @param Y, outcome vector
#' @param X, data.frame of variables that Y is a function of.
#' @param A, treatment vector
#' @param W, vector of variables A is a function of.
#' @param newdata, dataframe of X, stacked with X when A=1 and X when A=0, 
#' in that order
#' @param method, the SuperLearner meta learning method
#' @param SL.library, SuperLearner Library for finding outcome model
#' @param SL.libraryG, SuperLearner Library for the treatment mechanism
#' @param V, the number of folds
#' @param mc.cores,  number of cores to use for parallel processing the 
#' SuperLearner. Note, this parallelizes across the folds not within 
#' SuperLearner
#' @return A list with 5 elements:
#' initdata: the initdata argument for running tmle with gentmle function
#' 
#' Qcoef: the avg SuperLearner coef for each model in the outcome regression
#' 
#' Gcoef: the avg SuperLearner coef for each model in the treatment mech 
#' regression
#' 
#' Qrisk: the avg SuperLearner risk for each model in the outcome regression
#' 
#' Grisk: the avg SuperLearner risk for each model in the treatment mech
#' 
#' inds: the indices for all the val sets, stacked to match 
#' 
#' @export
#' @example /inst/examples/exampleATEandBV.R
SL.stack = function(Y, X, A, W, newdata, method, SL.library, 
                    SL.libraryG, V=10, mc.cores = 1, ...) {
  # 
  # X = X
  # Y = data$Y
  # A = data$A
  # W = X
  # W$A = NULL
  # newdata = newdata
  # method = "method.NNloglik"
  n = length(Y)
  folds = make_folds(n, V=V)
  stack = mclapply(folds, FUN = function(x) {
    # x=folds[[5]]
    tr = x$training_set
    val = x$validation_set
    nt=length(tr)
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
sim_hal = function(n, g0, Q0) {
  
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
  
  halresults <- hal(Y = data$Y,newX = newdata,
                    X = X, family = binomial(),
                    verbose = FALSE, parallel = FALSE)
  
  Qk = halresults$pred[1:n]
  Q1k = halresults$pred[n+1:n]
  Q0k = halresults$pred[2*n+1:n]
  
  initest = var(Q1k-Q0k)
  
  halresultsG <- hal(Y = data$A,newX = W,
                     X = W, family = binomial(),
                     verbose = FALSE, parallel = FALSE)
  gk = halresultsG$pred[1:n]
  Qcoef = Gcoef = 0
  
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
  
  
  steps = c(sigma_info$steps, sigmait_info$steps, simul_info$steps, simuljl_info$steps,
            simuljer_info$steps, simul_info$converge, ATE_info$steps)
  converge = c(sigma_info$converge, sigmait_info$converge,simul_info$converge, 
               simuljl_info$converge, simuljer_info$converge, simul_info$converge, 
               ATE_info$converge)
  
  ci_sig = ci_gentmle(sigma_info)[c(2,4,5)]
  ci_sigit = ci_gentmle(sigmait_info)[c(2,4,5)]
  ci_simul = ci_gentmle(simul_info)[2,c(2,4,5)]
  ci_simuljl = ci_gentmle(simuljl_info)[2,c(2,4,5)]
  ci_simuljer = ci_gentmle(simuljer_info)[2,c(2,4,5)]
  
  ci_simulATE = ci_gentmle(simul_info)[1,c(2,4,5)]
  ci_ATE = ci_gentmle(ATE_info)[c(2,4,5)]  
  
  cis = c(ci_sig,ci_sigit, ci_simul, ci_simuljl, ci_simuljer,ci_simulATE,ci_ATE)
  names(converge) = names(steps) = names(cis)[c(1,4,7,10,13,16,19)]=
    c("sig", "sigit", "simul", "simul_line", "simul_full","simulATE","ATE")
  results = c(cis, initest = initest,initest_ATE = initest_ATE, 
              steps = steps, converge = converge)
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
SL.stack1 = function(Y, X, A, W, newdata, method, SL.library, SL.libraryG, 
                     cv = TRUE, V=10, ...) {
  # 
  # X = X
  # Y = data$Y
  # A = data$A
  # W = X
  # W$A = NULL
  # newdata = newdata
  # method = "method.NNloglik"
  n = length(Y)
  if  (!cv) V = 1
  folds = make_folds(n, V=V)
  stack = lapply(folds, FUN = function(x) {
    # x=folds[[5]]
    if (!cv) {tr = val = 1:n} else {
    tr = x$training_set
    val = x$validation_set
  }
    nt=length(tr)
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
                Qrisk = Qrisk, Grisk = Grisk, inds = val))
  })
  
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
# 
 
#' @title sim_cv
#' @description Function that simulates data and performs TMLE estimates
#' data is 4 covariates and binary treatment and outcome.  The covariates
#' are generated according to the gendata function.
#' @param n, sample size
#' @param g0, treatment mechanism formula 
#' @param Q0, outcome model formula
#' @param SL.library, SuperLearner library for outcome predictions
#' @param SL.libraryG, SuperLearner Library for treatment mechanism
#' @param cv, set to TRUE for CV-TMLE
#' @param single, always set to FALSE.
#' @return  a vector with the following elements in this order:
#' TMLE confidence intervals each with estimate, left and right bounds
#' one step single parameter TMLE for blip variance, the iterative one
#' step TMLE for blip variance, simultaneous one step TMLE (which 
#' estimates average treatment effect and blip variance) for blip 
#' variance, simultaneous iterative TMLE "line" option
#' (estimates average treatment effect and blip variance) for blip 
#' variance, simultaneous iterative TMLE "full" option
#' (estimates average treatment effect and blip variance) for blip 
#' variance. All of the same CI's just mentioned except for initial
#' estimate that is main terms and interactions glm excep the "line"
#' and "full" models. one step simultaneous TMLE (estimating average
#' treatment effect and blip variance) for average treatment effect,
#' TMLE just for Average Treatment Effect, the same as the previous
#' two CI's where initial estimates for the TMLE's are glm with 
#' main terms and interactions. initial estimate for blip variance
#' using SuperLearner, initial estimate for blip variance using glm 
#' with main terms and interactions, initial estimate for ATE using
#' superlearner, initial estimate for ATE using glm with main terms
#' and interactions, steps to convergence, logical of whether the 
#' algorithm converged, Outcome SuperLearner coefficients, treatment
#' mechanism SuperLearner coefficients, Outcome SuperLearner risk,
#' Propensity score SuperLearner risk
#' @export
#' @example /inst/examples/example_sim_cv.R
sim_cv = function(n, g0, Q0, SL.library, SL.libraryG, method = "method.NNLS", 
                  cv = TRUE, single = FALSE) {
  data = gendata(n, g0, Q0)
  
  X = data
  X1 = X0 = X
  X0$A = 0
  X1$A = 1
  Y = data$Y
  A = data$A
  W = X
  W$A = NULL
  W$Y = NULL
  
  if (single) {} else{
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
  stack = SL.stack1(Y=Y, X=X, A=A, W=W, newdata=newdata, method=method, 
                      SL.library=SL.library, SL.libraryG=SL.libraryG,cv = cv)
  # proc.time() - time
  
  initdata = stack$initdata
  }
  
  # stack
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

# input data.frame with A, Y and covariates spit out lr CI based on delta method
# remember the order of vars is A, mainterms, then interactions
#' @export
LR.BVinference = function(n, g0, Q0) {
  
  data = gendata(n, g0 = g0, Q0 = Q0)  
  
  # set up treatment and outcome plus stack for pred
  data1 = data0 = data
  data1$A = 1
  data0$A = 0
  
  newdata = rbind(data,data1,data0)
  
  mainform = paste0(paste(colnames(data)[2:4],"+",collapse=""),colnames(data)[5])
  mainform

  main.int = paste0("Y~A*(",mainform,")")
  main.int = formula(main.int) 
  main.int
  
  newdata = model.matrix(main.int,newdata)
  newdata = as.data.frame(newdata[,-1])
  colnames(newdata)[2:ncol(newdata)] = paste0("X",2:ncol(newdata))
  
  # fit the regression
  X = newdata[1:n,]
  X$Y = data$Y
  Qfit = glm(Y~.,data=X,
             family='binomial')
  # predictions over data, A=1 and A=0
  Qk = predict(Qfit,type='response')
  Q1k = predict(Qfit,newdata=newdata[(n+1):(2*n),],type='response')
  Q0k = predict(Qfit,newdata=newdata[(2*n+1):(3*n),],type='response')
  
  # covariates and treatment for convenient use
  X = newdata
  X$Y = NULL
  X=cbind(int = rep(1,n),X)
  head(X)
  # calculate the score
  score_beta = sapply(1:n,FUN = function(x) {
    X[x,]*(data[x,"Y"]-Qk[x])
  })
  
  # averaging hessians to approx the deriv of hessian and mean then inverse
  hessian = lapply(1:n,FUN = function(x) {
    mat = -(1-Qk[x])*Qk[x]*as.numeric(X[x,])%*%t(as.numeric(X[x,]))
    return(mat)
  })
  fisher = -Reduce('+', hessian)/n
  M = solve(fisher)
  
  # calculate the IC for beta
  IC_beta = apply(score_beta,2,FUN = function(x) M%*%as.numeric(x))
  
  # SE_test = apply(IC_beta,1,sd)*sqrt(n-1)/n
  # SE_test
  
  blip = Q1k-Q0k
  ate = mean(Q1k-Q0k)
  # calculate the deriv to mult by IC_beta
  deriv1 = rowMeans(vapply(1:n, FUN = function(x) {
    return((1-Q1k[x])*Q1k[x]*as.numeric(X[(n+x),])-(1-Q0k[x])*Q0k[x]*
                              as.numeric(X[(2*n+x),]))
  }, FUN.VALUE=rep(1,10)))
  
  deriv = rowMeans(vapply(1:n, FUN = function(x) {
    return(2*(blip[x]-ate)*((1-Q1k[x])*Q1k[x]*as.numeric(X[(n+x),])-(1-Q0k[x])*Q0k[x]*
             as.numeric(X[(2*n+x),])))
  }, FUN.VALUE=rep(1,10)))
  
  psi = var(blip)
  # connect both parts of IC to form the full one
  
  IC = apply(IC_beta,2,FUN = function(x) t(deriv)%*%x) + (blip - ate)^2 - psi
  IC1 = apply(IC_beta, 2, FUN = function(x) t(deriv1)%*%x) + blip -ate
  # standard error
  SE = sd(IC)*sqrt((n-1))/n
  SE1 = sd(IC1)*sqrt((n-1))/n
  CI = c(psi_bv=psi,left=psi-1.96*SE,right=psi+1.96*SE)
  
  CI_ate = c(psi_ate = ate, ate-1.96*SE1,right=ate+1.96*SE1)
  
  corM = cor(data.frame(IC=IC, IC1=IC1))
  corM
  
  Z = rmvnorm(1000000,c(0,0),corM)
  zabs = apply(Z,1,FUN = function(x) max(abs(x)))
  zscore = quantile(zabs,.95)
  CI_simul_ate = c(psi_ate_simul = ate, left = ate - zscore*SE1, right = ate + zscore*SE1)  
  CI_simul_bv = c(psi_bv_simul = psi, left = psi - zscore*SE, right = psi + zscore*SE) 

  return(c(CI,CI_simul_bv,CI_ate,CI_simul_ate))
}

#' @export
noise.1 = function(data,n,rate, biasQ, sdQ) {
  rnorm(nrow(data), 
        with(data, biasQ(1,W1,W2,W3,W4,n=n,rate=rate)),
        with(data, sdQ(1,W1,W2,W3,W4,n=n,rate=rate)))
}

#' @export
# noise on barQ(0,W) is correlated per draw so not any more variant
noise.0 = function(data,noise,n,rate, biasQ, sdQ) {
  .5*noise + 
    sqrt(.75)*rnorm(nrow(data), 
                    with(data, biasQ(0,W1,W2,W3,W4,n=n,rate=rate)),
                    with(data, sdQ(0,W1,W2,W3,W4,n=n,rate=rate)))
}

#' @export
noise = function(data, noise_1, noise_0, biasQ, sdQ) {
  with(data, A*noise_1+(1-A)*noise_0)
}

#' @export
simBlipvar = function(n,rate, g0, Q0, biasQ, sdQ){
  # tack on the noise and use as an initial estimate
  data = gendata_noise(n, g0, Q0)
  noise_1 = noise.1(data, n, rate, biasQ, sdQ)
  noise_0 = noise.0(data, noise_1, n, rate, biasQ, sdQ)
  noise_A = noise(data, noise_1, noise_0 ,biasQ, sdQ)
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

#' @export
getRes = function(allresults,B, ATE0, var0, varind = c(4,10), 
                  ateind = c(7,11)) {
  # allresults = L[[3]]
  # B=1000
  results = vapply(1:B,FUN = function(x) (allresults[[x]]),
                   FUN.VALUE=allresults[[1]])
  results1=as.matrix(results)
  results=apply(results1,2,as.numeric)
  row.names(results)=row.names(results1)
  results=t(results)
  # results = results[results[,"converges.sigma"]==0,]
  cov.sig.it = mean(results[,2]<=var0&var0<=results[,3])
  cov.sig.1step = mean(results[,5]<=var0&var0<=results[,6])
  cov.ate = mean(results[,8]<=ATE0&ATE0<=results[,9])
  # cov.onestepEst = mean(results[,36]<=var0&var0<=results[,37])
  coverage = c(cov.sig.it=cov.sig.it,
               cov.sig.1step=cov.sig.1step,
               cov.ate=cov.ate)
  
  ests=c(results[,7],results[,11])
  type = c(rep("TMLE",B),rep("init",B))
  
  ateests = data.frame(ests=ests,type=type)
  ggover1 = ggplot(ateests,aes(ests, fill=type)) + 
    geom_density(alpha=.3)+
    scale_fill_manual(values=c("red", "blue"))+
    theme(axis.title.x = element_blank())+ggtitle("ATE sampling distributions")+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size=8))
  ggover1 = ggover1+geom_vline(xintercept = c(ATE0),color=c("green"))+
    # annotate("text",x=c(.2,.3),y=0,label=c("Event1","Event2"),hjust=1,angle=-40)
    # geom_text(x=.1,y=-5,label="asdfs",angle=-40,size = 3,hjust=1)
    geom_vline(xintercept=mean(results[,7]),color = "red")+
    geom_vline(xintercept=mean(results[,11]),color = "blue")+
    geom_vline(xintercept=ATE0,color = "green")
  # breaks = round(c(seq(min(ests),max(ests),length.out=4),ATE0,
  #            mean(results[,11]),mean(results[,7])),3)
  # ggover1 + scale_x_continuous(labels = c(round(seq(min(ests),max(ests),length.out=4),3),"ATE0",
  #                                         "mean init", "mean of TMLE")[order(breaks)],
  #                              breaks=breaks[order(breaks)])
  ggover1=ggdraw(add_sub(ggover1,"Truth at green line", 
                         x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                         vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                         colour = "black", size = 8, angle = 0, lineheight = 0.9))
  
  ests = c(results[,4],results[,10])
  type= c(rep("one step",B),rep("init",B))
  # ests = c(results[,23],results[,7],results[,34])
  # type= c(rep("LR",B),rep("one-step multi",B),rep("init",B))
  varests = data.frame(ests=ests,type=type)[1:B,]
  
  ggover2 = ggplot(varests,aes(ests)) + 
    geom_density(alpha=.3)+
    theme(axis.title.x = element_blank())+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size=8,angle=315))+
    ggtitle("blip variance sampling distributions")
  ggover2 = ggover2+geom_vline(xintercept = var0,color="green")+
    geom_vline(xintercept=mean(results[,4]),color = "black")+
    geom_vline(xintercept=mean(results[,10]),color = "red")
  # breaks = round(c(seq(min(varests[,1]),max(varests[,1]),length.out=4),mean(results[,10])),3)
  # ggover2 = ggover2 + scale_x_continuous(labels = c(round(seq(min(varests[,1]),max(varests[,1]),length.out=4),3),
  #                                         "mean init")[order(breaks)],
  #                              breaks=breaks[order(breaks)])
  ggover2=ggdraw(add_sub(ggover2,"truth at green line. mean of initial est at red line",x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                         vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                         colour = "black", size = 8, angle = 0, lineheight = 0.9))
  
  performance.sig = t(apply(results[,varind], 2, perf,var0))
  performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
  
  res = list(performance.ate=performance.ate, performance.sig=performance.sig, 
             coverage=coverage,ggover1=ggover1, ggover2=ggover2)
  return(res)
}

#' @export
noise_analysis = function(n, rate, truth, coverage, Q0, biasQ, sdQ) {
  # rate = .33
  # # V=1
  # # truth = gendata(10000)
  # n
  Q1_true = with(truth, Q0(1,W1,W2,W3,W4))
  Q0_true = with(truth, Q0(0,W1,W2,W3,W4))
  Q_true = with(truth,Q0(A,W1,W2,W3,W4))
  # g_true = with(truth, g0(W1,W2,W3,W4))
  
  noise.1 = function(data,n,rate) {
    rnorm(nrow(data), with(data, biasQ(1,W1,W2,W3,W4,n=n,rate=rate)), 
          with(data, sdQ(1,W1,W2,W3,W4,n=n,rate=rate)))}
  # noise on barQ(0,W) is correlated per draw so not any more variant
  noise.0 = function(data,noise,n,rate) {
    .5*noise +
      sqrt(.75)*rnorm(nrow(data), with(data, biasQ(0,W1,W2,W3,W4,n=n,rate=rate)),
                    with(data, sdQ(0,W1,W2,W3,W4,n=n,rate=rate)))}
  
  noise = function(data, noise_1, noise_0) with(data, A*noise_1+(1-A)*noise_0)
  
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
  
  # hist(Q_test)
  # hist(Q_true)
  # hist(blip_test)
  # hist(blip_true)
  # hist(g_test)
  ate_bias = mean(blip_test) - ATE0
  var_bias = var(blip_test) - var0
  
  title = paste("sample size ",n, "variance bias = ",round(var_bias,5))
  df = data.frame(blip = c(blip_test,blip_true), type = c(rep("test",1e6),rep("true",1e6)))
  gg_testvstrue = ggplot(df,aes(x=blip,color=type))+geom_density()+ggtitle(title)+
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size=8,angle=315)) 
  caption = paste0("coverage = ",coverage)
  gg_testvstrue = ggdraw(add_sub(gg_testvstrue,caption, x= .05, y = 0.5, hjust = 0, vjust = 0.5, vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                                 colour = "black", size = 9, angle = 0, lineheight = 0.9))
  # var_bias
  # ate_bias
  
  results = list(L2_blip = n^.25*L2_blip, L2_Q = n^.25*L2_Q,
                 ate_bias = ate_bias, var_bias = var_bias,plot = gg_testvstrue)
  return(results)
}
