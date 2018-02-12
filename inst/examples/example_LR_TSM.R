
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4) + W3 +W4)

# get a randomly drawn dataframe under the specified model
data = gendata(1000, g0_linear, Q0_linear)
big = gendata(1e6, g0_linear, Q0_linear)
# get the truth setting A to 1
setA = 1
truth = mean(with(big, Q0_linear(A=setA,W1,W2,W3,W4)))
truth
# well-specified model
Qform = formula("Y ~ W1 + W2 + A*(W3 + W4)")

# specifying the covariates, treatment and outcome
W = data[,2:5]
A = data$A
Y = data$Y

undebug(LR.TSM)
undebug(IC.beta)
# should cover each truth 95 percent of the time.
info = LR.TSM(W=W,A=A,Y=Y,Qform=Qform, setA = 1, alpha = .05)

# get CI
info$CI
# get influence curve
info$IC


