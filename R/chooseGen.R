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
  # n = 1000; d = 1; pos = .05; minATE = .05; minBV = .05; depth = 1; maxterms = 1; minterms = 1; mininters = 1
  # num.binaries = floor(d/4)
  
  if (minterms == 0) 
    stop("minterms must be atleast 1")
  if (mininters > minterms) 
    stop("minimum number of interactions cannot exceed number of covariate terms, obviously, hello!!!")
  
  # sample size of population
  N = 1e+06
  
  # randomly drawn binary distributions 
  poss.binaries = sample(1:d, num.binaries)
  r = runif(num.binaries, 0.3, 0.7)
  
  # the basic matrix of confounders
  Wmat = vapply(1:d, FUN = function(col) {
    if (col <= num.binaries) {
      return(rbinom(N, 1, r[col]))
    } else {
      return(rnorm(N, 0, 1))
    }
  }, FUN.VALUE = rep(1,N))
  
  # All of the interaction combos
  choos = lapply(1:depth, FUN = function(x) {
    c = combn(1:d, x)
    if (!is.matrix(c)) c = as.matrix(c)
    return(c)
  })
  
  # select the interaction terms to be included according to maxterms and minterms
  s = 0
  if (d == 1) {
    dfG = as.matrix(types[[sample(1:4, 1)]](Wmat))
  } else {
    while (s <= minterms) {
      terms = lapply(choos, FUN = function(x) {
        no.terms = sample(0:min(maxterms, ncol(x)),1)
        select.cols = sample(1:ncol(x), no.terms)
        return(select.cols)
      })
      s = sum(unlist(lapply(terms, FUN = function(x) length(x))))
    }
    
    # interact the columns and store in list 
    col.comb = lapply(1:length(terms), FUN = function(a) {
      col.choos = terms[[a]]
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][, x]
        return(rowProds(Wmat, cols = col.inds))
      }, FUN.VALUE = rep(1, N))
    })
    
    # put the cols in one matrix
    dfG = do.call(cbind, col.comb)
    
    # types of transformations to apply, can add many more
    types = list(function(x) sin(x), function(x) cos(x), 
                 function(x) x^2, function(x) x)
    
    # transform the columns by plugging into randomly drawn functions and standardize
    # so no variable dominates unnecessarily
    dfG = apply(dfG, 2, FUN = function(col) {
      if (all(col == 1 | col ==0)) {
        v = (col - mean(col))/sd(col)
        return(v)
      } else {
        v = types[[sample(1:4, 1)]](col)
        v = (v - mean(v))/sd(v)
        return(v)
      }
    })
  }
  
  skewage = runif(1, -1, 1)  
  dfG = cbind(dfG, rep(skewage, N))
  coef_G = runif(ncol(dfG), -1, 1)
  # satisfying positivity constraints
  tol = TRUE
  its = 0
  while (tol & its < 20) {
    PG0 = plogis(dfG %*% coef_G)
    coef_G = 0.8 * coef_G
    its = its + 1
    tol = mean(PG0 < pos) > 0.01 | mean(PG0 > 1 - pos) > 
      0.01
  }
  
  # Creating A and standardizing
  PG0 = gentmle2::truncate(PG0, pos)
  # hist(PG0, breaks = 100)
  A = rbinom(N, 1, PG0)
  
  # for barQ first select interaction terms for just the W's
  s = 0
  if (d == 1) {
    onlyW = types[[sample(1:4, 1)]](Wmat)
    dfQW = matrix((onlyW - mean(onlyW))/sd(onlyW), ncol = 1)
    dfQWA = cbind(dfQW, (A - mean(A)/sd(A)))
    if (mininters == 1) {
      df_inter = dfQW
      fcn = types[[sample(1:4, 1)]]
      df_interA = fcn(dfQW*A)
      means = mean(df_interA)
      sds = sd(df_interA)
      df_interA = (df_interA - means)/sds
      df_inter = (fcn(dfQW) - means)/sds
      df_inter0 = (fcn(rep(0,N)) - means)/sds
      dfQ = cbind(dfQWA, df_interA)
      dfQ1 = dfQWA
      dfQ0 = dfQWA
      dfQ1[,2] = (1 - mean(A))/sd(A)
      dfQ0[,2] = -mean(A)/sd(A)
      dfQ1 = cbind(dfQ1, df_inter)
      dfQ0 = cbind(dfQ0, df_inter0)
    }
  } else {
    while (s <= minterms) {
      termsQW = lapply(choos, FUN = function(x) {
        no.terms = sample(0:min(maxterms, ncol(x)),1)
        select.cols = sample(1:ncol(x), no.terms)
        return(select.cols)
      })
      s = sum(unlist(lapply(terms, sum)))
    }
    
    # for barQ select interaction terms of W's that will interact with A
    s = 0
    while (s <= mininters) {
      terms_inter = lapply(choos, FUN = function(x) {
        no.terms = sample(0:min(maxterms, ncol(x)),1)
        select.cols = sample(1:ncol(x), no.terms)
        return(select.cols)
      })
      s = sum(unlist(lapply(terms, sum)))
    }
    
    # interact the columns and store in list 
    col.combQ = lapply(1:length(termsQW), FUN = function(a) {
      col.choos = termsQW[[a]]
      dfQW = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][, x]
        return(rowProds(Wmat, cols = col.inds))
      }, FUN.VALUE = rep(1, N))
    })
    
    # put the cols in one matrix
    dfQWA = do.call(cbind, col.combQ)
    dfQWA = cbind(dfQWA, A)
    
    col.comb_inter = lapply(1:length(terms_inter), FUN = function(a) {
      col.choos = terms_inter[[a]]
      dfQ_inter = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][, x]
        return(rowProds(Wmat, cols = col.inds))
      }, FUN.VALUE = rep(1, N))
    })
    
    # put the cols in one matrix
    dfQ_inter = do.call(cbind, col.comb_inter)
    dfQ_interA = apply(dfQ_inter, 2, FUN = function(col) A*col)
    
    dfQWA = apply(dfQWA, 2, FUN = function(col) {
      if (all(col == 1 | col ==0)) {
        v = (col - mean(col))/sd(col)
        return(v)
      } else {
        v = types[[sample(1:4, 1)]](col)
        v = (v - mean(v))/sd(v)
        return(v)
      }
    })
    
    ftypes = lapply(1:ncol(dfQ_interA), FUN = function(col) {
      if (all(dfQ_interA[,col] == 1 | dfQ_interA[,col] ==0)) {
        return(types[[4]])
      } else {
        v = types[[sample(1:4, 1)]]
        return(v)
      }
    })
    
    dfQ_interA = vapply(1:length(ftypes), FUN = function(col) ftypes[[col]](dfQ_interA[,col]), FUN.VALUE = rep(1,N))
    dfQ_inter = vapply(1:length(ftypes), FUN = function(col) ftypes[[col]](dfQ_inter[,col]), FUN.VALUE = rep(1,N))
    dfQ_inter0 = vapply(1:length(ftypes), FUN = function(col) ftypes[[col]](rep(0,N)), FUN.VALUE = rep(1,N))
    
    means = apply(dfQ_interA, 2, FUN = function(col) mean(col))
    sds = apply(dfQ_interA, 2, FUN = function(col) sd(col))
    
    dfQ_interA = apply(dfQ_interA, 2, FUN = function(col) (col - mean(col))/sd(col))
    dfQ_inter = vapply(1:ncol(dfQ_inter), FUN = function(col) {
      (dfQ_inter[,col] - means[col])/sds[col]
    }, FUN.VALUE = rep(1,N))
    dfQ_inter0 = vapply(1:ncol(dfQ_inter0), FUN = function(col) {
      (dfQ_inter0[,col] - means[col])/sds[col]
    }, FUN.VALUE = rep(1,N))
    
    dfQ = cbind(dfQWA, dfQ_interA)
    dfQW1 = dfQWA
    dfQW1[, ncol(dfQW1)] = (1 - mean(A))/sd(A)
    dfQ1 = cbind(dfQW1, dfQ_inter)
    
    dfQW0 = dfQWA
    dfQW0[, ncol(dfQW0)] = - mean(A)/sd(A)
    dfQ0 = cbind(dfQW0, dfQ_inter0)
  }
  
  TXpos = ncol(dfQWA)  
  a = runif(1, 0, 1)
  coef_Q = runif(ncol(dfQ), -a, a)
  PQ1 = plogis(dfQ1 %*% coef_Q)
  PQ0 = plogis(dfQ0 %*% coef_Q)
  blip_true = PQ1 - PQ0
  ATE0 = mean(blip_true)
  BV0 = var(blip_true)

  jj = 1
  while (abs(ATE0) <= minATE & jj <= 20) {
    coef_Q[TXpos] = 1.2*coef_Q[TXpos]
    PQ1 = plogis(dfQ1 %*% coef_Q)
    PQ0 = plogis(dfQ0 %*% coef_Q)
    blip_true = PQ1 - PQ0
    ATE0 = mean(blip_true)
    jj = jj + 1
  } 
  
  jj = 1 
  if (ncol(dfQ_inter) != 0) {
    while (BV0 <= minBV & jj <= 20) {
      coef_Q[(ncol(dfQWA) + 1):ncol(dfQ)] = 1.2 * coef_Q[(ncol(dfQWA) + 1):ncol(dfQ)]
      PQ1 = plogis(dfQ1 %*% coef_Q)
      PQ0 = plogis(dfQ0 %*% coef_Q)
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
  Wn = as.data.frame(Wmat[S, ])
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
