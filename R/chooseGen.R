# We have 4 covariates and 1 treatment, which can be binary or continuous, 
# we first choose which are binary and which are continuous (rnorms)

# declare the dimension of W and n

# now we have W1, W2, W3, W4 which are random variables as plain binomials and fcns of normals
# we can now interact these fcns and create new variables

#' @title get.dgp
#' @description randomly creates a dgp, for now only 4 dimensional of confounders, binary 
#' treatment.
#' @param n, sample size
#' @param d, set to 4 only
#' @param pos, a small value to make sure prop scores are between in (pos, 1 - pos)
#' @param minBV, minimum blip variance
#' @param depth, specify depth of interaction--must be less than d
#' @param maxterms, maximum terms per interaction.  For example, this would limit 
#' two way interactions to maximally 10 terms as well as three way or main terms.
#' With high dimension it is wise to set this low because it might take a while
#' otherwise.  Still in development, be cautious.
#' @param minterms sets a minimum number of total covariate terms, including 
#' interactions with eachother--do not set lower than 1.
#' @param mininters sets the minimum number of interactions with treatment to include
#' This must be bigger or equal to minterms
#' @return  a sample DF, the true average treatment effect, ATE0 and blip variance
#' BV0, the sample pscores, PGn, the sample true blips, blip_n, the sample 
#' true prob of death under treatment, PQ1n, and prob of death under control
#' PQ0n
#' @export
#' @example /inst/examples/example_get.dgp.R
get.dgp = function(n, d, pos = .01, minBV = 0, depth, maxterms, minterms, mininters) {
  # n=1000; d=8;pos=.01;minBV=.03; depth = 4; maxterms = 10; minterms = d; mininters = 2
  # sample size for getting the truth
  
  if (minterms == 0) stop("minterms must be atleast 1")
  
  if (mininters > minterms) stop("minimum number of interactions cannot exceed number of covariate terms, obviously, hello!!!")
  
  N = 1e6
  
  # choose binaries--at most 1 out of 4 covariates is a binary
  num.binaries = floor(d/4)
  poss.binaries = sample(1:d, num.binaries)
  binaries = rep(FALSE, d)
  binaries[poss.binaries] = c(as.logical(rbinom(num.binaries, 1, .5)))
  binaries = (1:d)[binaries]
  binaries
  
  # randomly select a binary prob between .3 and .7
  r = runif(1, .3, .7)
  # randomly select chisq, normals, betas, uniforms--normalized
  
  Wbig = matrix(rep(NA, d*N), ncol = d)
  for (a in 1:d) {
    if (a %in% binaries) {
      V = rbinom(N, 1, r)
    } else {
      V = rnorm(N, 0, 1)
    }
    Wbig[,a] = V 
  }
  
  Wmat = Wbig
  Wbig = as.data.frame(Wbig)
  colnames(Wbig) = paste0("W", 1:d)
  
  # get contins cols
  conts = vapply(1:d, FUN = function(x) !(x %in% binaries), FUN.VALUE = TRUE)
  contins  = (1:d)[conts]
  
  #####
  #####
  # we use the previous parameters to generate the transformed vars and truth
  #####
  #####
  choos = lapply(1:depth, FUN = function(x) {
    combn(1:d, x)
  })
  # create the transformed covariates--coeffs to be added
  dfs = lapply(1:2, FUN = function(j){
    # vary functional forms
    # j = 1
    types = list(function(x) sin(x), function(x) cos(x), function(x) x^2, 
                 function(x) x)
    
    # choosing from combos of 2, 3 and whether we have 4 way interaction
    # This should be generalized to larger dimensions but for now only 4
    s = 0
    while (s <= minterms) {
      terms = lapply(choos, FUN = function(x) {
        no.terms = sample(0:min(maxterms, ncol(x)))
        select.cols = sample(1:ncol(x), no.terms)
        return(1:ncol(x) %in% select.cols)
      })
      s = sum(unlist(lapply(terms, sum)))
    }
    
    df = vapply(1:d, FUN = function(i) {
      if (i %in% binaries) {
        return(types[[4]](Wmat[,i]))
      } else {
        return(types[[sample(1:4, 1)]](Wmat[,i]))
      }
    }, FUN.VALUE = rep(1, N))
    
    df_final = df[,terms[[1]]]
    
    for (a in 2:length(terms)) {
     if (sum(terms[[a]] != 0)) {
       col.choos = which(terms[[a]])
       df = vapply(col.choos, FUN = function(x) {
         col.inds = choos[[a]][,x]
         return(rowProds(Wmat, cols = col.inds))
     }, FUN.VALUE = rep(1, N)) 
       df_final = cbind(df_final, df)
     }
    }
    
    df_final = apply(df_final, 2, FUN = function(x) {
      (x - mean(x))/sd(x)
    })
    return(df_final) 
  })
  
  
  skewer = sample(1:10, 1)
  if (skewer <= 4) p = -.1 else {if (skewer <= 7) p = 0 else p = .1}
  # p
  df_final = dfs[[1]] + p
  # for (a in 1:ncol(dfs[[1]])) hist(dfs[[1]][,a], breaks = 200)
  # generating coeffs for prop score 
  skewfactor = sample(2:5,1)
  a = 1
  coef_G = runif(ncol(df_final), -a, skewfactor*a)
  # coef_G
  # assure no true practical positivity violations without having neverending loop
  # pscores = plogis(dfs[[1]] %*% coef_G)
  tol = TRUE
  its = 0
  while (tol & its < 20) {
    PG0  = plogis(df_final %*% coef_G)
    coef_G = .8*coef_G
    its = its + 1
    tol = mean(PG0 < pos) > .01 | mean(PG0 > 1-pos) > .01
  }
  # hist(PG0)
  # max(PG0)
  # min(PG0)
  # its
  # skewfactor
  # p
  PG0 = gentmle2::truncate(PG0, pos)
  # hist(PG0, breaks = 100)
  
  # now get the A for everyone
  A = rbinom(N, 1, PG0)
  
  # setting up dataframe for both Q0k found previously and Bk, the interactions terms
  
  no.inters = sample(mininters:(ncol(dfs[[2]])+1), 1)
  select.cols = sample(1:(ncol(dfs[[2]])+1), no.inters)
  inters = 1:(ncol(dfs[[2]])) %in% select.cols
  
  
  dfinter0 = vapply(which(inters), FUN = function(col) dfs[[2]][,col]*A, FUN.VALUE = rep(1, N))
  
  if ((ncol(dfs[[2]])+1) %in% select.cols) {
    TX = (A - mean(A))/sd(A)
    dfTX = cbind(dfs[[2]], TX)
    TXcol = dfTX[,ncol(dfTX)]
    TX0 = -mean(A)/sd(A)
    TX1 = (1 - mean(A))/sd(A)
    dfQ0 = cbind(dfs[[2]], TX0)
    dfQ1 = cbind(dfs[[2]], TX1)
    dfQ = cbind(dfTX, dfinter0)
    dfQ1 = cbind(dfQ1, dfs[[2]][, which(inters)])
  } else {
    dfTX = dfQ0 = dfQ1 = dfs[[2]]
    dfQ = cbind(dfTX, dfinter0)
    dfQ1 = cbind(dfQ1, dfs[[2]][, which(inters)])
  }



  # choosing coeffs for the Q generator
  a = 1
  coef_Q = runif(ncol(dfQ), -a, a)
  
  # The following helps get the blip var over .025
  BV0 = 0
  
  jj = 1
  # coef_Q
  while (BV0 <= minBV & jj < 20) {
    coef_Q[(ncol(dfTX) + 1):ncol(dfQ)] = 1.2*coef_Q[(ncol(dfTX) + 1):ncol(dfQ)]
    PQ1  = plogis(dfQ1 %*% coef_Q)
    PQ0 = plogis(dfQ0 %*% coef_Q[1:ncol(dfQ0)])
    blip_true = PQ1 - PQ0
    ATE0 = mean(blip_true)
    BV0 = var(blip_true)
    
    jj = jj + 1
  }
  # coef_Q
  # BV0
  # now that we have our coeffs, get true probs of sample and draw Y
  PQ  = plogis(dfQ %*% coef_Q)
  # hist(PQ0, breaks = 200)
  # BV0
  Y = rbinom(N, 1, PQ)
  
  # mean(Y*A/PG0)-mean(Y*(1-A)/(1-PG0))
  # FF = data.frame(Y=Y,A=A,PG0=PG0)
  # mean(with(FF, Y*A/PG0-Y*(1-A)/(1-PG0)))

  # now we select 1000 out of this population for W, A, Y
  S = sample(1:N, n)
  PQ1n  = PQ1[S]
  PQ0n = PQ0[S]
  PGn = PG0[S]
  blip_n = PQ1n - PQ0n
  An = A[S]
  Yn = Y[S]
  Wn = Wbig[S,]
  DF = cbind(Wn, An, Yn)
  colnames(DF)[c((d+1),(d+2))] = c("A", "Y")
  
  # return the sample true blips and sample barQ1, barQ0, as well as truths and the DF
  return(list(BV0 = BV0, ATE0 = ATE0, DF = DF, blip_n = blip_n, PQ1n = PQ1n, PQ0n = PQ0n, PGn = PGn))
  
}

# testDF = get.dgp(n = 1000, d = 12, pos = .01, minBV = .03, depth = 4, maxterms = d, minterms = d, mininters = 12)
# 
# head(chubs$DF)
# chubs$BV0
# chubs$ATE0
# hist(chubs$blip_n,100)
# hist(chubs$PQ1n,100)
# hist(chubs$PQ0n,100)
# hist(chubs$PGn,100)
# max(chubs$PGn)
# min(chubs$PGn)

# big = gendata(1e6, g0_1, Q0_2)
# gtrue = with(big, g0_1(W1,W2,W3,W4))
# gtrue[1:10]
# mean(big$A*big$Y/gtrue - (1 - big$A)*big$Y/(1 - gtrue))
# get.truth(g0_1, Q0_2)
