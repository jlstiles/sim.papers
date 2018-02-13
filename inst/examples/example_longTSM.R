###
# The following example is for pt treatment
###
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

data = gendata(n, g0_linear, Q0_linear)

formulas = list(formula("Y ~ W1 + W2 + A*(W3 + W4)"))
Ynodes = c("Y")
Anodes = c("A")
setA = 1

TSMinfo = long.TSM(data = data, Ynodes = Ynodes, Anodes = Anodes, 
                   formulas = formulas, setA = setA, alpha = .05)

# get CI
TSMinfo$CI
# get influence curve
TSMinfo$IC

####
# example with longitudinal intervention
####

data(longdata)
head(data)

Ynodes = c("Y_1", "Y_2","Y_3","Y_4")
Anodes = c("A1_0", "A1_1","A1_2","A1_3")

# specify the formulas
formula0 = formula("Y_1 ~ L2_0 + L1_0 + A1_0") 
formula1 = formula("Y_2 ~ L2_1 + A1_1")
formula2 = formula("Y_3 ~ L2_2 + A1_2")
formula3 = formula("Y_4 ~ L2_3 + A1_3")
formulas = list(formula0, formula1, formula2, formula3)
setA = c(0,1,1,1)

TSMinfo = long.TSM(data = data, Ynodes = Ynodes, Anodes = Anodes, 
                   formulas = formulas, setA = setA)

TSMinfo$CI
TSMinfo$IC
