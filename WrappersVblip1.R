###################
###################
# gam
###################
###################

# Three screeners
# choose main terms
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

screen10 = function (Y, X, family, method = "pearson", rank = 10, ...) 
{
  listp <- vapply(colnames(X), FUN = function(x) {
    value = ifelse(var(X[,x]) <= 0, 1, cor.test(X[,x], y = Y, method = method)$p.value)
    if (x == "A") value = 0
    return(value)
  },FUN.VALUE = 9)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}

environment(screen10) <- asNamespace("SuperLearner")

screen.Main = function (Y, X, family, method = "pearson", rank = 6, ...) 
{
  whichVariable <- rep(FALSE,ncol(X))
  whichVariable[1:5] <- rep(TRUE,5)
  return(whichVariable)
}

environment(screen.Main) <- asNamespace("SuperLearner")

SL.gam3 <- function (Y, X, newX, family, obsWeights, deg.gam = 3, cts.num = 10,...)
{
  SuperLearner:::.SL.require("gam")
  if ("mgcv" %in% loadedNamespaces())
    warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~", paste(paste("s(",
                                                    colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,
                                                    ")", sep = ""), collapse = "+"), "+", paste(colnames(X[,
                                                                                                           !cts.x, drop = FALSE]), collapse = "+")))
  }
  else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(",
                                                    colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,
                                                    ")", sep = ""), collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X),
                                              collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = X, family = family,
                      control = gam::gam.control(maxit = 50, bf.maxit = 50),
                      weights = obsWeights)
  pred <- gam::predict.gam(fit.gam, newdata = newX, type = "response")
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam")
  return(out)
}

environment(SL.gam3) <- asNamespace("SuperLearner")

###########################
###########################
# XGB
###########################
###########################

create.Learner("SL.xgboost",tune=list(max_depth=1,minobspernode=3,shrinkage=.001,ntrees=10000))
xgbFull = SL.xgboost_1
create.Learner("SL.xgboost",tune=list(max_depth=4,minobspernode=3,shrinkage=.001,ntrees=2500))
xgbMain = SL.xgboost_1
create.Learner("SL.xgboost",tune=list(max_depth=2,minobspernode=6,shrinkage=.005,ntrees=1000))
xgb6 = SL.xgboost_1
create.Learner("SL.xgboost",tune=list(max_depth=1,minobspernode=6,shrinkage=.01,ntrees=2500))
xgb10 = SL.xgboost_1

#######################
#######################
# nnet
#######################
#######################
# for mainterms, size 5, decay=.1, full model: size =2, decay .1

nnetMain = function (Y, X, newX, family, obsWeights, size = 5, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .1)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}

environment(nnetMain) <- asNamespace("SuperLearner")

nnetMain.2 = function (Y, X, newX, family, obsWeights, size = 2, ...)
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

environment(nnetMain.2) <- asNamespace("SuperLearner")

nnetMain.3 = function (Y, X, newX, family, obsWeights, size = 4, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE,
                           maxit = 500, linout = FALSE, weights = obsWeights,decay = .3)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}

environment(nnetMain.3) <- asNamespace("SuperLearner")


nnetFull = function (Y, X, newX, family, obsWeights, size = 2, ...)
{
  SuperLearner:::.SL.require("nnet")
  if (family$family == "gaussian") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE,
                           trace = FALSE, maxit = 500, weights = obsWeights)
  }
  if (family$family == "binomial") {
    fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE,
                           maxit = 500, linout = FALSE, weights = obsWeights, decay = .1)
  }
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}

environment(nnetFull) <- asNamespace("SuperLearner")

#######################
#######################
# rpart
#######################
#######################
rpartPrune = function (Y, X, newX, family, obsWeights, cp = 0.001, minsplit = 5,
          xval = 5, maxdepth = 15, minbucket = 5, ...)
{
  SuperLearner:::.SL.require("rpart")
  if (family$family == "gaussian") {
    fit.rpart <- rpart(Y ~ ., data = data.frame(Y, X), control = rpart.control(cp = cp,
                                                                               minsplit = minsplit, xval = xval, maxdepth = maxdepth,
                                                                               minbucket = minbucket), method = "anova", weights = obsWeights)
    CP <- fit.rpart$cptable[which.min(fit.rpart$cptable[,
                                                        "xerror"]), "CP"]
    fitPrune <- prune(fit.rpart, cp = CP)
    pred <- predict(fitPrune, newdata = newX)
  }
  if (family$family == "binomial") {
    fit.rpart <- rpart(Y ~ ., data = data.frame(Y, X), control = rpart.control(cp = cp,
                                                                               minsplit = minsplit, xval = xval, maxdepth = maxdepth,
                                                                               minbucket = minbucket), method = "class", weights = obsWeights)
    CP <- fit.rpart$cptable[which.min(fit.rpart$cptable[,
                                                        "xerror"]), "CP"]
    fitPrune <- prune(fit.rpart, cp = CP)
    pred <- predict(fitPrune, newdata = newX)[, 2]
  }
  fit <- list(object = fitPrune, fit = fit.rpart, cp = CP)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.rpart")
  return(out)
}

environment(rpartPrune) <- asNamespace("SuperLearner")

#######################
#######################
# ranger
#######################
#######################
# for choose 10 nodesize 10 and 1000 trees with mtry 3 up to 2500 for full interaction,
# 1000 for choose 6

create.Learner("SL.ranger",tune=list(min.node.size = 10, mtry = 3, num.trees=2500))
rangerFull = SL.ranger_1
create.Learner("SL.ranger",tune=list(min.node.size = 10, mtry = 3, num.trees=2500))
ranger10 = SL.ranger_1

#######################
#######################
# glmnet
#######################
#######################

create.Learner("SL.glmnet",tune=list(alpha = c(0,.5,1)))

#######################
#######################
# earth
#######################
#######################

earthMain = function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
                        nk = max(21, 2 * ncol(X) + 1), ...)
{
  SuperLearner:::.SL.require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree,
                              nk = nk, penalty = penalty)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree,
                              nk = nk, penalty = penalty,
                              minspan=10,
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}
environment(earthMain) <- asNamespace("SuperLearner")

earthFull = function (Y, X, newX, family, obsWeights, id, degree = 1, penalty = 3,
                      nk = max(21, 2 * ncol(X) + 1), ...)
{
  SuperLearner:::.SL.require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree,
                              nk = nk, penalty = penalty)
  }
  if (family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree,
                              nk = nk, penalty = penalty,
                              minspan=10,
                              glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}
environment(earthFull) <- asNamespace("SuperLearner")

#######################
#######################
# glm
#######################
#######################
glm.mainint = function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  mainform = paste0(paste(colnames(X)[2:4],"+",collapse=""),colnames(X)[5])
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

#######################
#######################
# stepAIC
#######################
#######################
