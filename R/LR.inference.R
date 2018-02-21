#' @title IC.beta
#' @description computes IC for beta coefficients of logistic regression
#' mainly a helper function.  
#' @param W, matrix or data.frame of covariates
#' @param A, a binary vector of treatment assignments
#' @param Y, a binary vector of outcomes
#' @param Qform, a formula for Y in terms of the covariates as input in glm
#' 
#' @return  a list with elements IC_beta and Qfit, the glm fit object.  
#' @export
IC.beta = function(data,OC=NULL, Ynode, Anodes, Qform, verbose = FALSE, parallelize = FALSE) {
  n = nrow(data)
  # This option is for when feeding in sequential regression
  if (!is.null(OC)) data[,Ynode] = OC
  # we only fit on non deaths or uncensored
  cens = is.na(data[,Ynode])
  data = data[!cens,]
  n1 = nrow(data)
  # form the design matrix based on the formula
  X = model.matrix(Qform,data)
  X = as.data.frame(X[,-1])
  # fit the regression
  Y = data[,Ynode]
  if (!verbose) {
    Qfit = suppressWarnings(stats::glm(Y~.,data=X,
                    family='binomial'))
  } else {
    Qfit = stats::glm(Y~.,data=X,family='binomial')
  }
  goods = 1:(ncol(X)+1)
  if (any(is.na(coef(Qfit)))) {
    print(paste0("you have a singular covariance matrix so we will refit without these variables",
                 paste(names(coef(Qfit))[is.na(coef(Qfit))], collapse = " ")))
    goods = which(!is.na(coef(Qfit)))
    X = X[,(goods-1)]
    if (!verbose) {
      Qfit = suppressWarnings(stats::glm(Y~.,data=X,
                                         family='binomial'))
    } else {
      Qfit = stats::glm(Y~.,data=X,family='binomial')
    }
  }
  # predictions over data
  Qk = predict(Qfit,type='response')
  
  X$Y = NULL
  X=cbind(int = rep(1,n1),X)
  
  # calculate the score
  # score_beta = sapply(1:n1,FUN = function(x) {
  #   X[x,]*(Y[x]-Qk[x])
  # })
  
  if (parallelize) cores = getOption("mc.cores",parallel::detectCores()) else cores = 1L
  score_beta = mclapply(1:n1,FUN = function(x) {
    X[x,]*(Y[x]-Qk[x])
  }, mc.cores = getOption("mc.cores", cores))
  
  LL = length(score_beta[[1]])
  score_beta = vapply(1:length(score_beta), FUN = function(x){
    return(unlist(score_beta[[x]]))
           }, FUN.VALUE = rep(1,LL))
  
  # averaging hessians to approx the true average then invert as per IC
  hessian = mclapply(1:n1,FUN = function(x) {
    mat = -(1-Qk[x])*Qk[x]*as.numeric(X[x,])%*%t(as.numeric(X[x,]))
    return(mat)
  }, mc.cores = getOption("mc.cores", cores))
  
  # M1 = summary(Qfit)$cov.unscaled*n1
  fisher = -Reduce('+', hessian)/n1
  M = solve(fisher)
  
  Xfull = matrix(rep(NA,(n*ncol(X))),nrow = n)
  Xfull = as.data.frame(Xfull)
  Xfull[!cens,] = X
  colnames(Xfull) = colnames(X)
  # calculate the IC for beta
  IC_beta = matrix(rep(0, nrow(M)*n), nrow = nrow(M))
  IC_beta[,!cens] = apply(score_beta,2,FUN = function(x) M%*%as.numeric(x))
  IC_beta = IC_beta
  return(list(IC_beta = IC_beta, Qfit = Qfit, X = Xfull, hessian = M, goods = goods))
  
}

#' @title long.TSM
#' @description computes IC and psi for logistic regression plug-in estimator of treatment
#' specific mean for survival data in wide form. IN DEVELOPMENT
#' @param data, data.frame of variables in time ordering from left to right
#' @param Ynodes, character vector of time-ordered Ynodes
#' @param Anodes, character vector of time-ordered Ynodes
#' @param formulas, list of formulas for the conditional means
#' @param setA the value to which you intervene on A,vector of length that of Anodes
#' @param alpha significance level for two-sided CI.
#' @return  a list with elements, CI for the confidence interval and  IC for the 
#' influence curve.
#' @export
#' @example /inst/examples/example_longTSM.R
long.TSM = function(data, Ynodes, Anodes, formulas, setA, alpha = .05, parallelize = FALSE)
{
  n = nrow(data)
  Yinds = vapply(Ynodes, FUN = function(x) grep(x,colnames(data)), FUN.VALUE = 1)
  times = 1:length(setA)
  times = times[order(times,decreasing = TRUE)]
  for (t in times){
    # If t is the end time we just do a regression on the outcome
    if (t == max(times)) {
      IC_tplus1 = 0
      OC = NULL
      if (t == 1) design = data else design = data[,-Yinds[1:(t-1)]]
      ICinfo_t = IC.beta(data = design, OC = OC, Ynode = Ynodes[t], 
                         Anode = Anodes[1:t], Qform = formulas[[t]], parallelize = parallelize)
      IC_t = ICinfo_t$IC_beta
      # otherwise we recursively proceed to define the IC
    } else {
      IC_tplus1 = IC_t
      ICinfo_tplus1 = ICinfo_t
      # get the indice to make the design matrix based on the formula
      Yind = grep(Ynodes[t+1], colnames(data))
      # for the outcomes at t which are 0 we use prev regression o.w use Y = 1
      # goods are those that did not die
      Y_t = data[,Ynodes[t]]
      goods = vapply(Y_t, FUN = function(x) {
        t = ifelse(!is.na(x), x==0, FALSE) 
      }, FUN.VALUE = TRUE)
      reals = vapply(Y_t, FUN = function(x) {
        ifelse(!is.na(x), x==1, FALSE) 
      }, FUN.VALUE = TRUE)
      # make the outcome predictions based on previous beta by grabbing previous 
      # design, intervening then setting up the design--should prob 
      # switch to datatable
      Xa_tplus1 = data[goods,-Yinds[1:t]]
      for (i in 1:(t+1)) {
        col = grep(Anodes[i], colnames(Xa_tplus1))
        Xa_tplus1[,col] = setA[i]
      }
      Xa_tplus1 = model.matrix(formulas[[t+1]],Xa_tplus1)
      Xa_tplus1 = Xa_tplus1[,(ICinfo_tplus1$goods)]
      OC = rep(NA,n)
      OC[goods] = plogis(Xa_tplus1 %*% ICinfo_tplus1$Qfit$coef)
      OC[reals] = 1
      # get the new beta
      # if (t == 1) design = data else design = data[,-Yinds[1:(t-1)]]
      ICinfo_t = IC.beta(data = data, OC = OC, Ynode = Ynodes[t], 
                         Anode = Anodes[1:t], Qform = formulas[[t]],parallelize = parallelize)
      X_t = ICinfo_t$X[goods,]
      OC = OC[goods]
      # This is to create the M matrix from the paper
      hess = mclapply(1:length(OC),FUN = function(x) {
        mat = (1-OC[x])*OC[x]*as.numeric(X_t[x,])%*%t(as.numeric(Xa_tplus1[x,]))
        return(mat)
      }, mc.cores = getOption("mc.cores", parallel::detectCores()))
      M = Reduce('+', hess)/length(OC)
      M = ICinfo_t$hessian %*% M 
      # form the IC
      IC_temp = apply(IC_tplus1,2,FUN = function(x) M%*%as.numeric(x))
      IC_t = ICinfo_t$IC_beta + IC_temp
      
    }
  } 
  
  n = nrow(data)
  # if (t == 1) XA = data else XA = data[,-Yinds[1:(t-1)]]
  XA = data
  XA[,Anodes[1]] = setA[1]
  XA = model.matrix(formulas[[1]],XA)
  XA = XA[,(ICinfo_t$goods)]
  QAk = plogis(XA %*% ICinfo_t$Qfit$coef)
  
  score = rowMeans(sapply(1:n,FUN = function(x) {
    QAk[x]*(1 - QAk[x])*XA[x,]
  })) 
  
  psi = mean(QAk)
  IC = apply(IC_t, 2, FUN = function(x) sum(score*x)) + QAk - psi
  SE = sd(IC)*sqrt(n-1)/n
  
  qq = qnorm(1-alpha/2)
  CI = c(psi = psi, left = psi - qq*SE, right = psi + qq*SE)
  
  return(list(CI = CI, IC = IC))
}



#' @title LR.inference
#' @description Function that gives inference for logistic regression plug-in
#' estimators of ATE and Blip Variance.  
#' @param W, matrix or data.frame of covariates
#' @param A, a binary vector of treatment assignments
#' @param Y, a binary vector of outcomes
#' @param Qform, a formula for Y in terms of the covariates as input in glm
#' @param alpha, significance level for the (1-alpha)100 percent CI's. 0.05 is default
#' @param simultaneous.inference, TRUE if user wants simultaneous confidence
#' bounds for both ATE and blip variance at level alpha. default is FALSE
#' 
#' @return  if simultaneous.inference is specified as TRUE then will return a vector giving
#' pt estimate, left and right bound for ATE, simultaneous ATE CI, blip variance, 
#' and simultaneous blip variance.  Otherwise gives pt estimate, left and right bound for ATE
#' and blip variance.  
#' @export
#' @example /inst/examples/example_LR_inference.R
LR.inference = function(W, A, Y, Qform, alpha = .05, simultaneous.inference = FALSE) {
  n = length(Y)
  
  X = as.data.frame(cbind(A,W,Y))
  X0 = X1 = X
  X0$A = 0
  X1$A = 1
  
  newdata = rbind(X, X1, X0)
  newdata = model.matrix(Qform,newdata)
  newdata = as.data.frame(newdata[,-1])
  colnames(newdata)[2:ncol(newdata)] = paste0("X",2:ncol(newdata))
  
  # fit the regression
  Qfit = stats::glm(Y~.,data=newdata[1:n,],
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
    X[x,]*(Y[x]-Qk[x])
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
  }, FUN.VALUE=rep(1,ncol(X))))
  
  deriv = rowMeans(vapply(1:n, FUN = function(x) {
    return(2*(blip[x]-ate)*((1-Q1k[x])*Q1k[x]*as.numeric(X[(n+x),])-(1-Q0k[x])*Q0k[x]*
                              as.numeric(X[(2*n+x),])))
  }, FUN.VALUE=rep(1,ncol(X))))
  
  psi = var(blip)
  # connect both parts of IC to form the full one
  
  IC = apply(IC_beta,2,FUN = function(x) t(deriv)%*%x) + (blip - ate)^2 - psi
  IC1 = apply(IC_beta, 2, FUN = function(x) t(deriv1)%*%x) + blip -ate
  # standard error
  SE = sd(IC)*sqrt((n-1))/n
  SE1 = sd(IC1)*sqrt((n-1))/n
  
  qq = qnorm(1-alpha/2)
  CI = c(bv_delta = psi, left = psi - qq*SE, right = psi + qq*SE)
  
  CI_ate = c(ate_delta = ate, left = ate - qq*SE1, right = ate + qq*SE1)
  
  if (simultaneous.inference) {
    corM = stats::cor(data.frame(IC=IC, IC1=IC1))
    Z = rmvnorm(1000000,c(0,0),corM)
    zabs = apply(Z,1,FUN = function(x) max(abs(x)))
    zscore = quantile(zabs, 1-alpha)
    CI_simul_ate = c(ate_deltasimul = ate, left = ate - zscore*SE1, right = ate + zscore*SE1)  
    CI_simul_bv = c(bv_deltasimul = psi, left = psi - zscore*SE, right = psi + zscore*SE) 
    
    return(c(CI_ate, CI_simul_ate, CI, CI_simul_bv))
  } else {
    return(c(CI_ate, CI))
  }
}

#' @ export
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
