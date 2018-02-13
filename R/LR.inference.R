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
IC.beta = function(data,OC=NULL, Ynode, Anodes, Qform) {
  n = nrow(data)
  if (!is.null(OC)) data[,Ynode] = OC
  cens = is.na(data[,Ynode])
  data = data[!cens,]
  n1 = nrow(data)
  X = model.matrix(Qform,data)
  X = as.data.frame(X[,-1])
  colnames(X) = paste0("X",1:(ncol(X)))
  # fit the regression
  Y = data[,Ynode]
  Qfit = stats::glm(Y~.,data=X,
                    family='binomial')
  # predictions over data, A=1 and A=0
  Qk = predict(Qfit,type='response')
  
  X$Y = NULL
  X=cbind(int = rep(1,n1),X)
  # head(X)
  # calculate the score
  score_beta = sapply(1:n1,FUN = function(x) {
    X[x,]*(Y[x]-Qk[x])
  })
  
  # averaging hessians to approx the deriv of hessian and mean then inverse
  hessian = lapply(1:n1,FUN = function(x) {
    mat = -(1-Qk[x])*Qk[x]*as.numeric(X[x,])%*%t(as.numeric(X[x,]))
    return(mat)
  })
  fisher = -Reduce('+', hessian)/n1
  M = solve(fisher)
  
  # calculate the IC for beta
  IC_beta = matrix(rep(0, nrow(M)*n), nrow = nrow(M))
  IC_beta[,!cens] = apply(score_beta,2,FUN = function(x) M%*%as.numeric(x))
  IC_beta = IC_beta*n/n1
  return(list(IC_beta = IC_beta, Qfit = Qfit, X = X, hessian = M))
  
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
long.TSM = function(data, Ynodes, Anodes, formulas, setA, alpha = .05)
{
  times = 1:length(setA)
  times = times[order(times,decreasing = TRUE)]
  for (t in times){
    # undebug(IC.beta)
    if (t == max(times)) {
      IC_tplus1 = 0
      OC = NULL
      ICinfo_t = IC.beta(data = data, OC = OC, Ynode = Ynodes[t], 
                         Anode = Anodes[1:t], Qform = formulas[[t]])
      IC_t = ICinfo_t$IC_beta
    } else {
      IC_tplus1 = IC_t
      ICinfo_tplus1 = ICinfo_t
      n0 = nrow(ICinfo_tplus1$X)
      Yind = grep(Ynodes[t+1], colnames(data))
      Y_t = data[,Ynodes[t]]
      goods = vapply(Y_t, FUN = function(x) {
        t = ifelse(!is.na(x), x==0, FALSE) 
      }, FUN.VALUE = TRUE)
      
      Xa_tplus1 = data[goods,1:Yind]
      for (i in 1:(t+1)) {
        col = grep(Anodes[i], colnames(Xa_tplus1))
        Xa_tplus1[,col] = setA[i]
      }
      Xa_tplus1 = model.matrix(formulas[[t+1]],Xa_tplus1)
      OC = rep(NA,n)
      OC[goods] = plogis(Xa_tplus1 %*% ICinfo_tplus1$Qfit$coef)
      ICinfo_t = IC.beta(data = data, OC = OC, Ynode = Ynodes[t], 
                         Anode = Anodes[1:t], Qform = formulas[[t]])
      X_t = ICinfo_t$X
      OC = OC[!is.na(OC)]
      hess = lapply(1:n0,FUN = function(x) {
        mat = (1-OC[x])*OC[x]*as.numeric(X_t[x,])%*%t(as.numeric(Xa_tplus1[x,]))
        return(mat)
      })
      M = Reduce('+', hess)/n0
      M = ICinfo_t$hessian %*% M 
      IC_temp = apply(IC_tplus1,2,FUN = function(x) M%*%as.numeric(x))
      IC_t = ICinfo_t$IC_beta + IC_temp
      
    }
  } 
  
  
  n = nrow(data)
  XA = data
  XA[,Anodes[1]] = setA[1]
  XA = model.matrix(formulas[[1]],XA)
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
