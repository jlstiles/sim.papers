devtools::install_github("jlstiles/Simulations", force = TRUE)
# devtools::install_github("jlstiles/halplus", force = TRUE)
# devtools::install_github("benkeser/halplus", force = TRUE)

library(Simulations)
get.truth(g0_linear, Q0_trig1)

set.seed(14819)
n=200
data =gendata(n, g0=g0_linear,Q0=Q0_trig1)
X=data
X$Y=NULL
Y=data$Y
X0=X1=X
X0$A=0
X1$A=1
newdata = rbind(X,X1,X0)
# 
# undebug(SL.hal)
Qfit = SuperLearner(Y=Y, X=X, newX = newdata, family = binomial(), method = "method.NNloglik",
             SL.library = c("SL.hal","SL.glm"), cvControl = list(V=10L))

Qfit$Z
Qfit$cvControl
Qfit$library.predict
Qfit$metaOptimizer
Qfit$cvRisk

# 
# cl = makeCluster(detectCores(), type = "SOCK")
# registerDoSNOW(cl)
time = proc.time()
halres <- hal(Y = Y, newX = newdata, X = X, family =  binomial(),
                                   verbose = FALSE, parallel = FALSE)

halres1 <- hal(Y = Y, newX = newdata, X = X, family = gaussian(),
              verbose = FALSE, parallel = FALSE)
timehal = proc.time() - time

Q = halres$pred[1:n]
Q1 = halres$pred[n+1:n]
Q0 = halres$pred[2*n+1:n]
risk = mean(Y*log(Q) + (1-Y)*log(1-Q))
risk
max(Q)
min(Q)

Q_1 = halres1$pred[1:n]
Q1_1 = halres1$pred[n+1:n]
Q1_0 = halres1$pred[2*n+1:n]


esthal = var(Q1-Q0)
esthal1 = var(Q1_1-Q1_0)

esthalATE = mean(Q1-Q0)
esthal1ATE = mean(Q1_1-Q1_0)

save(esthal, esthal1, esthalATE, esthal1ATE, risk, file = "test.RData")

# time = proc.time()
# halres9001 <- fit_hal(Y = Y, X = X, family = 'binomial')
# timehal9001 = proc.time() - time
# 
# QkH = predict(halres9001, new_data = newdata, type = 'response')[1:n]
# Q1kH = predict(halres9001, new_data = newdata)[n+1:n]
# Q0kH = predict(halres9001, new_data = newdata)[2*n+1:n]
# 
# riskhal9001 = mean(Y*log(QkH)+(1-Y)*log(1-QkH))
# 
# esthal9001 = var(Q1kH-Q0kH)
# esthal9001ATE = mean(Q1kH-Q0kH)
# 
# timehal
# timehal9001
# 
# # answers are different
# esthal
# esthal9001
# 
# esthalATE
# esthal9001ATE
# 
# riskhal
# riskhal9001
# bit::bit
