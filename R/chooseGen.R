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
get.dgp = function(n, d, pos = 0.01, minATE = -2, minBV = 0, depth, maxterms, minterms, 
                    mininters, num.binaries = floor(d/4)) 
{
  # n = 1000; d = 4; pos = .05; minBV = .05; depth = 4; maxterms = 2; minterms = 2; mininters = 2
  # num.binaries = floor(d/4)
  if (minterms == 0) 
    stop("minterms must be atleast 1")
  if (mininters > minterms) 
    stop("minimum number of interactions cannot exceed number of covariate terms, obviously, hello!!!")
  
  # sample size of population
  N = 1e+06
  
  # make at most 1/4 of vars binaries for now. 
  poss.binaries = sample(1:d, num.binaries)
  r = runif(num.binaries, 0.3, 0.7)

  Wmat = vapply(1:d, FUN = function(col) {
    if (col <= num.binaries) {
      return(rbinom(N, 1, r[col]))
    } else {
      return(rnorm(N, 0, 1))
    }
  }, FUN.VALUE = rep(1,N))
  
  choos = lapply(1:depth, FUN = function(x) {
    c = combn(1:d, x)
    if (!is.matrix(c)) c = as.matrix(c)
    return(c)
  })
  
  
  dfs = lapply(1:2, FUN = function(j) {
    types = list(function(x) sin(x), function(x) cos(x), 
                 function(x) x^2, function(x) x)
    s = 0
    if (d == 1) {
      df_final = apply(Wmat, 2, types[[sample(1:4, 1)]])
    } else {
      while (s <= minterms) {
        terms = lapply(choos, FUN = function(x) {
          no.terms = sample(0:min(maxterms, ncol(x)),1)
          select.cols = sample(1:ncol(x), no.terms)
          return(1:ncol(x) %in% select.cols)
        })
        s = sum(unlist(lapply(terms, sum)))
      }}
      
    col.comb = lapply(1:length(terms), FUN = function(a) {
      col.choos = which(terms[[a]])
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][, x]
        return(rowProds(Wmat, cols = col.inds))
      }, FUN.VALUE = rep(1, N))
    })
    
    df_final = do.call(cbind, col.comb)
    
    df_final = apply(df_final, 2, FUN = function(col) {
      if (all(col == 1 | col ==0)) {
        v = (col - mean(col))/sd(col)
        return(v)
      } else {
        v = types[[sample(1:4, 1)]](col)
        v = (v - mean(v))/sd(v)
        return(v)
        }
      })
    
    return(df_final)
  })
  skewer = sample(1:10, 1)
  if (skewer <= 4){ 
    p = -0.1
  } else {
    if (skewer <= 7) 
      p = 0
    else p = 0.1
  }
  df_final = dfs[[1]] + p
  skewfactor = sample(2:5, 1)
  a = 1
  coef_G = runif(ncol(df_final), -a, skewfactor * a)
  tol = TRUE
  its = 0
  while (tol & its < 20) {
    PG0 = plogis(df_final %*% coef_G)
    coef_G = 0.8 * coef_G
    its = its + 1
    tol = mean(PG0 < pos) > 0.01 | mean(PG0 > 1 - pos) > 
      0.01
  }
  PG0 = gentmle2::truncate(PG0, pos)
  A = rbinom(N, 1, PG0)
  TX = (A - mean(A))/sd(A) 
  
  no.inters = sample(mininters:ncol(dfs[[2]]), 1)
  select.cols = sample(1:ncol(dfs[[2]]), no.inters)
  dfinter0 = vapply(select.cols, FUN = function(col) dfs[[2]][,col] * A, FUN.VALUE = rep(1, N))
  
  dfQ0 = cbind(dfs[[2]], rep(-mean(A)/sd(A), N))
  dfQ = cbind(dfs[[2]], TX, dfinter0)
  dfQ1 = cbind(dfs[[2]], rep(1-mean(A)/sd(A), N), dfs[[2]][, select.cols])
  
  a = runif(1, 0, 1)
  coef_Q = runif(ncol(dfQ), -a, a)
  BV0 = 0
  jj = 1
  while (ATE0 <= minBV & jj < 20) {
    coef_Q[ncol(dfQ)] = 1.2*coef_Q[ncol(dfQ)]
    coef_Q[(ncol(dfQ0) + 1):ncol(dfQ)] = 1.2 * coef_Q[(ncol(dfQ0) + 1):ncol(dfQ)]
    PQ1 = plogis(dfQ1 %*% coef_Q)
    PQ0 = plogis(dfQ0 %*% coef_Q[1:ncol(dfQ0)])
    blip_true = PQ1 - PQ0
    jj = jj + 1
  } 
  
  jj = 0 
  if (no.inters != 0) {
    while (BV0 <= minBV & jj < 20) {
      coef_Q[(ncol(dfQ0) + 1):ncol(dfQ)] = 1.2 * coef_Q[(ncol(dfQ0) + 1):ncol(dfQ)]
      PQ1 = plogis(dfQ1 %*% coef_Q)
      PQ0 = plogis(dfQ0 %*% coef_Q[1:ncol(dfQ0)])
      blip_true = PQ1 - PQ0
      BV0 = var(blip_true)
      jj = jj + 1
    } 
  } else {
    PQ1 = plogis(dfQ1 %*% coef_Q)
    PQ0 = plogis(dfQ0 %*% coef_Q[1:ncol(dfQ0)])
    blip_true = PQ1 - PQ0
    BV0 = var(blip_true)
  }
  
  PQ = plogis(dfQ %*% coef_Q)
  Y = rbinom(N, 1, PQ)
  Y[PQ <= .00001] = 0
  Y[PQ >= .99999] = 1
  S = sample(1:N, n)
  PQ1n = PQ1[S]
  PQ0n = PQ0[S]
  PQn = PQ[S]
  PGn = PG0[S]
  blip_n = PQ1n - PQ0n
  An = A[S]
  Yn = Y[S]
  Wn = Wmat[S, ]
  DF = cbind(Wn, An, Yn)
  colnames(DF)[c((d + 1), (d + 2))] = c("A", "Y")
  colnames(DF)[1:d] = paste0("W",1:d)
  return(list(BV0 = BV0, ATE0 = ATE0, DF = DF, blip_n = blip_n, 
              PQ1n = PQ1n, PQ0n = PQ0n, PQn = PQn, PGn = PGn))
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
