# trainxgb$results[order(trainxgb$results[,"RMSE"]),]
create.Learner("SL.xgboost", tune = list(shrinkage = .001, max_depth = 2, 
                                         minobspernode = 9, ntrees = 5000))
xgboostMain = SL.xgboost_1

create.Learner("SL.xgboost", tune = list(shrinkage = .001, max_depth = 1, 
                                         minobspernode = 9, ntrees = 5000))

xgboostFull = SL.xgboost_1

create.Learner("SL.xgboost", tune = list(shrinkage = .005, max_depth = 1, 
                                         minobspernode = 9, ntrees = 2500))
xgboostG = SL.xgboost_1

xg = create.Learner("SL.xgboost", tune = list(shrinkage = .01, max_depth = c(2,3), 
                                         minobspernode = 6, ntrees = 1000))
xgboost_2dG = SL.xgboost_1
xgboost_2d = SL.xgboost_2
# trainearth$results
# trainearthG$results
earth = create.Learner("SL.earth",tune=list(degree = c(1,2)))
earthFull = SL.earth_1
earth_2d = SL.earth_2

earth = create.Learner("SL.earth",tune=list(degree = 1))
earthFull = SL.earth_1
# use same for G
# and use stock SL.earth and SL.gam for the main terms model

create.Learner("SL.gam", tune = list(deg.gam = 1))
gamFull = SL.gam_1
# use Full for G as well as SL.gam

# tunennet$performances[order(tunennet$performances[,"error"]),]
# tunennetG$performances[order(tunennetG$performances[,"error"]),]

nnetMainG = function (Y, X, newX, family, obsWeights, size = 2, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .50)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}
environment(nnetMainG) <- asNamespace("SuperLearner")

nnetMainG1 = function (Y, X, newX, family, obsWeights, size = 2, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .2)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}
environment(nnetMainG1) <- asNamespace("SuperLearner")


nnetMain = function (Y, X, newX, family, obsWeights, size = 4, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE, skip = TRUE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .30)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}
environment(nnetMain) <- asNamespace("SuperLearner")

nnetMain1 = function (Y, X, newX, family, obsWeights, size = 2, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE, skip = TRUE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .05)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}
environment(nnetMain1) <- asNamespace("SuperLearner")

nnet_d2 = function (Y, X, newX, family, obsWeights, size = 4, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE, skip = TRUE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .01)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}
environment(nnetMain) <- asNamespace("SuperLearner")


# tunerf$performances[order(tunerf$performances[,"error"]),]
create.Learner("SL.ranger", tune = list(num.trees = 2500, mtry = 3, 
                                        min.node.size = 10))
rangerMain = SL.ranger_1

# tunerpart$performances[order(tunerpart$performances[,"error"]),]
create.Learner("SL.rpartPrune", tune = list(cp = .005, minsplit = 15,
                                            minbucket = 30))
rpartMain = SL.rpartPrune_1

screen6 = function (Y, X, family, method = "pearson", rank = 6, ...) 
{
  listp <- vapply(colnames(X), FUN = function(x) {
    value = ifelse(var(X[,x]) <= 0, 1, cor.test(X[,x], y = Y, method = method)$p.value)
    if (x == "A") value = 0
    return(value)
  },FUN.VALUE = 9)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}

environment(screen6) <- asNamespace("SuperLearner")

screen12 = function (Y, X, family, method = "pearson", rank = 12, ...) 
{
  listp <- vapply(colnames(X), FUN = function(x) {
    value = ifelse(var(X[,x]) <= 0, 1, cor.test(X[,x], y = Y, method = method)$p.value)
    if (x == "A") value = 0
    return(value)
  },FUN.VALUE = 9)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}

environment(screen12) <- asNamespace("SuperLearner")

screen.Main = function (Y, X, family, method = "pearson", rank = 8, ...) 
{
  whichVariable <- rep(FALSE,ncol(X))
  whichVariable[1:8] <- rep(TRUE,8)
  return(whichVariable)
}

environment(screen.Main) <- asNamespace("SuperLearner")

glm.mainint = function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  mainform = paste0(paste(colnames(X)[2:ncol(X)],"",collapse="+"))
  form = formula(paste0("Y ~", paste0(colnames(X)[1],"*(",mainform,")")))
  
  fit.glm <- glm(form, data = X, family = family, weights = obsWeights,
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}
environment(glm.mainint) <- asNamespace("SuperLearner")
