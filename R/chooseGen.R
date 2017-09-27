# We have 4 covariates and 1 treatment, which can be binary or continuous, 
# we first choose which are binary and which are continuous (rnorms)

# declare the dimension of W and n

# now we have W1, W2, W3, W4 which are random variables as plain binomials and fcns of normals
# we can now interact these fcns and create new variables

#' @export
get.info = function(n, d, pos = .01, minBV = .03) {
  # n=1000; d=4;pos=.01;minBV=.03
  # sample size for getting the truth
  N = 1e6
  # choose binaries--possibly 0 or 1 for now
  binaries = c(as.logical(rbinom(1, 1, .5)), rep(FALSE, 3))
  binaries = (1:d)[binaries]
  
  # randomly select a binary prob
  r = runif(1, .3, .7)
  
  Wbig = matrix(rep(NA, d*N), ncol = d)
  for (a in 1:d) {
    if (a %in% binaries) {
      V = rbinom(n, 1, r)
    } else {
      V = rnorm(n, 0, 1)
    }
    Wbig[,a] = V 
  }
  
  
  Wbig = as.data.frame(Wbig)
  colnames(Wbig) = paste0("W", 1:d)
  
  # get contins cols
  conts = vapply(1:d, FUN = function(x) !(x %in% binaries), FUN.VALUE = TRUE)
  contins  = (1:d)[conts]
  
  # vary functional forms
  types = list(function(x) sin(x), function(x) cos(x), function(x) x^2, 
               function(x) x)
  
  # get forms for both Q and G
  pars = list()
  for (F in 1:2) {
    # amplify the binary or not
    bin_coef = runif(1, -5, 5)
    
    # choosing functions to apply to each variable
    funclist = list()
    for (x in 1:d) {
      if (x %in% binaries) {
        funclist[[x]] = function(g) bin_coef*g
      } else {
        funclist[[x]] = types[[sample(1:4, 1)]]
      }
    }
    
    # choosing from combos of 2, 2 and whether we have 4 way interaction
    # This should be generalized to larger dimensions but for now only 4
    s = 0
    while (s <= 1) {
      choo2 = as.logical(rbinom(6, 1, .5))
      choo3 = as.logical(rbinom(4, 1, .5))
      way4 = rbinom(1, 1, .5)
      # choosing which main terms to include
      MTnum = sample(1:d, 1)
      MT = sample(1:d, MTnum)
      MT = MT[order(MT)]
      s = sum(choo2) + sum(choo3) + way4 + MTnum
    }
    
    # storing for later use
    pars[[F]] = list(choo2 = choo2, choo3 = choo3, way4 = way4, MT = MT, funclist = funclist, bin_coef = bin_coef)
  } 
  
  #####
  #####
  # we use the previous parameters to generate the transformed vars and truth
  #####
  #####
  
  # create the transformed covariates--coeffs to be added
  dfs = lapply(pars, FUN = function(x){
    N = 1e6
    df = vapply(1:d, FUN = function(i) {
      bin_coef = x$bin_coef
      funclist[[i]](Wbig[,i])
    }, FUN.VALUE = rep(1, N))
    
    choo2 = x$choo2
    choo3 = x$choo3
    way4 = x$way4
    MT = x$MT
    funclist = x$funclist
    
    # make number of interactions and which ones
    ways2 = matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), byrow = TRUE, nrow = 6)
    df_final = df[,MT]
    
    if (sum(choo2) != 0) {
      df_2way = vapply(which(choo2), FUN = function(combo) {
        df[, ways2[combo,1]]*df[, ways2[combo, 2]]
      }, FUN.VALUE = rep(1,N))
      df_final = cbind(df_final, df_2way)
    }
    
    # make number of 3 ways interactions
    ways3 = matrix(c(1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4), byrow = TRUE, nrow = 4)
    if (sum(choo3) != 0) {
      df_3way = vapply(which(choo3), FUN = function(combo) {
        df[, ways3[combo,1]]*df[, ways3[combo, 2]]*df[, ways3[combo, 3]]
      }, FUN.VALUE = rep(1,N))
      df_final = cbind(df_final, df_3way)
    }
    
    if (way4) {
      df_4way = df[, 1]*df[, 2]*df[, 3]*df[, 4]
      df_final = cbind(df_final, df_4way)
    }
    
    # 
    # if (all(c(all(!choo2), all(!choo3), !way4))) {
    #   df_final = df
    #   df_final1 = df1_true
    # }
    return(df_final) 
  })
  
  # generating coeffs for prop score 
  a = .7 
  coef_G = runif(ncol(dfs[[1]]), -a, a)
  
  # assure no true practical positivity violations without having neverending loop
  PG0  = plogis(dfs[[1]] %*% coef_G)
  maxG = max(PG0)
  minG = min(PG0)
  g_iter = 0
  while ((maxG >= (1 - pos) | minG <= pos) & g_iter < 10) {
    coef_G = .8*coef_G
    PG0  = plogis(dfs[[1]] %*% coef_G)
    maxG = max(PG0)
    minG = min(PG0)
    g_iter = g_iter+1
  }
  
  # now get the A for everyone
  A = rbinom(N, 1, PG0)
  
  # setting up dataframe for both Q0k found previously and Bk, the interactions terms
  s = 0
  while (s == 0) {
    inters = rbinom(ncol(dfs[[2]]), 1, .5)
    s = sum(inters)
  }
  
  dfinter0 = vapply(which(inters == 1), FUN = function(col) dfs[[2]][,col]*A, FUN.VALUE = rep(1, N))
  dfQ = cbind(dfs[[2]], dfinter0)
  dfQ1 = cbind(dfs[[2]], dfs[[2]][, inters])
  
  # choosing coeffs for the Q generator
  a = .7
  coef_Q = runif(ncol(dfQ), -a, a)
  
  # The following helps get the blip var over .025
  BV0 = 0
  C = 1
  mm = 5
  jj = 1
  while (BV0 <= minBV & jj < 11) {
    coef_Q[(ncol(dfs[[2]]) + 1):ncol(dfQ)] = C*coef_Q[(ncol(dfs[[2]]) + 1):ncol(dfQ)]
    PQ1  = plogis(dfQ1 %*% coef_Q)
    PQ0 = plogis(dfs[[2]] %*% coef_Q[1:ncol(dfs[[2]])])
    blip_true = PQ1 - PQ0
    ATE0 = mean(blip_true)
    BV0 = var(blip_true)
    C = 1.1*C
    jj = jj + 1
  }
  
  BV0
  # now that we have our coeffs, get true probs of sample and draw Y
  PQ  = plogis(dfQ %*% coef_Q)
  # hist(PQ_true0, breaks = 200)
  # BV0
  Y = rbinom(N, 1, PQ)
  
  # now we select 1000 out of this population for W, A, Y
  S = sample(1:N, n)
  PQ1n  = PQ1[S]
  PQ0n = PQ0[S]
  PGn = PG0[S]
  blip_n = PQ1n - PQ0n
  An = A[S]
  Yn = Y[S]
  Wn = Wbig[S,]
  DF = cbind(An, Wn, Yn)
  colnames(DF)[c(1,(d+2))] = c("A", "Y")
  
  # return the sample true blips and sample barQ1, barQ0, as well as truths and the DF
  return(list(BV0 = BV0, ATE0 = ATE0, DF = DF, blip_n = blip_n, PQ1n = PQ1n, PGn = PGn))
  
}

# info = get.info(n = 1000, d = 4)
# info$BV0
# info$ATE0
# mean(info$DF$A)
# mean(info$DF$Y)
# min(info$PGn)
# max(info$PGn)
# head(info$DF)
# hist(info$blip_n)
