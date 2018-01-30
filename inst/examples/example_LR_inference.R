# using built-in package functions, g0_linear and Q0_1 to specify
# pscore and outcome model probabilities
data = gendata(1000, g0_linear, Q0_1)

Qform = formula("Y ~ W2 + W3 + A*W4")

# specifying the covariates, treatment and outcome
W = data[,2:5]
A = data$A
Y = data$Y

info = LR.inference(W=W,A=A,Y=Y,Qform=Qform, alpha = .05)
info


