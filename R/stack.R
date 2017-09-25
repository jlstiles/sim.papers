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
SL.stack1 = function(Y, X, A, W, newdata, method, SL.library, SL.libraryG, 
                     cv = TRUE, V=10, SL = 10L, ...) {
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
    if (!cv) {
      tr = val = 1:n
    } else {
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
                      SL.library=SL.library, method = method,
                      id = NULL, verbose = FALSE, control = list(),
                      cvControl = list(V=SL), obsWeights = NULL)
    
    A = A[tr]
    W1 = W[tr,]
    newW = W[val,]
    gfit = SuperLearner(Y=A,X=W1,newX = newW, family = binomial(),
                        SL.library=SL.libraryG, method = "method.NNloglik", 
                        id = NULL, verbose = FALSE, control = list(),
                        cvControl = list(V=SL), obsWeights = NULL)
    
    
    if (length(gfit$coef[gfit$coef!=0])==1) {
      gk = gfit$library.predict[1:nv,gfit$coef!=0]
    } else {
      gk = gfit$library.predict[1:nv,gfit$coef!=0] %*% gfit$coef[gfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1) {
      Qk = Qfit$library.predict[1:nv,Qfit$coef!=0]
    } else {
      Qk = Qfit$library.predict[1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1) {
      Q1k = Qfit$library.predict[nv+1:nv,Qfit$coef!=0]
    } else {
      Q1k = Qfit$library.predict[nv+1:nv,Qfit$coef!=0] %*% Qfit$coef[Qfit$coef!=0]
    }
    
    if (length(Qfit$coef[Qfit$coef!=0])==1) {
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
  
  if (is.vector(Qcoef_mat)) {
    Qcoef = mean(Qcoef_mat)
    Qrisk = mean(Qrisk_mat)
  } else {
    Qcoef = rowMeans(Qcoef_mat)
    Qrisk = rowMeans(Qrisk_mat)
  }
  
  if (is.vector(Gcoef_mat)) {
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