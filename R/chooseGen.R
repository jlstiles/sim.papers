# We have 4 covariates and 1 treatment, which can be binary or continuous, 
# we first choose which are binary and which are continuous (rnorms)

# declare the dimension of W and n

# now we have W1, W2, W3, W4 which are random variables as plain binomials and fcns of normals
# we can now interact these fcns and create new variables

#' @export
get.info = function(n, d, truth) {

  # n=1000
  # d=4
  # truth=TRUE
  N = 1e6
  # choose binaries
  binaries = c(as.logical(rbinom(1, 1, .5)), rep(FALSE, 3))
  binaries = (1:d)[binaries]
  
  # make all variables, n  copies of covariates
  W = matrix(rep(NA, d*n), ncol = d)
  for (a in 1:d) {
    nombre = paste0("W", a)
    if (a %in% binaries) {
      r = runif(1, .3, .7)
      V = rbinom(n, 1, r)
    } else {
      V = rnorm(n, 0, 1)
    }
    W[,a] = V 
  }
  
  if (truth) {
    Wbig = matrix(rep(NA, d*N), ncol = d)
    for (a in 1:d) {
      nombre = paste0("W", a)
      if (a %in% binaries) {
        V = rbinom(n, 1, r)
      } else {
        V = rnorm(n, 0, 1)
      }
      Wbig[,a] = V 
    }
  }
  
  # We then choose functions for the continuous variables
  # We choose from trig functions, squares, 1st deg, categorical breaks
  
  # get conins cols
  conts = vapply(1:d, FUN = function(x) !(x %in% binaries), FUN.VALUE = TRUE)
  contins  = (1:d)[conts]

  types = list(function(x) sin(x), function(x) cos(x), function(x) x^2, 
               function(x) x)
  
  # get forms for both Q and G
  pars = list()
  for (F in 1:2) {
  # amplify the binary or not
  bin_coef = runif(1, -5, 5)
  
  funclist = list()
  for (x in 1:d) {
    if (x %in% binaries) {
      funclist[[x]] = function(g) bin_coef*g
    } else {
      funclist[[x]] = types[[sample(1:4, 1)]]
    }
  }
  
  choo2 = as.logical(rbinom(6, 1, .5))
  choo3 = as.logical(rbinom(4, 1, .5))
  way4 = rbinom(1, 1, .5)
  MTnum = sample(1:d, 1)
  MT = sample(1:d, MTnum)
  MT = MT[order(MT)]
  
  pars[[F]] = list(choo2 = choo2, choo3 = choo3, way4 = way4, MT = MT, funclist = funclist, bin_coef = bin_coef)
  } 
  
  dfs = lapply(pars, FUN = function(x){
    
    # make df of main terms and form fcns
    df = vapply(1:d, FUN = function(i) x$funclist[[i]](W[,i]), FUN.VALUE = rep(1, n))
    choo2 = x$choo2
    choo3 = x$choo3
    way4 = x$way4
    MT = x$MT
    funclist = x$funclist
    bin_coef = x$bin_coef
    # make number of interactions and which ones
    ways2 = matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), byrow = TRUE, nrow = 6)
    
    if (sum(choo2) != 0) {
      df_2way = vapply(which(choo2), FUN = function(combo) {
        df[, ways2[combo,1]]*df[, ways2[combo, 2]]
      }, FUN.VALUE = rep(1,n))
      df_final = cbind(df[,MT], df_2way)
    }
    
    # make number of 3 ways interactions
    ways3 = matrix(c(1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4), byrow = TRUE, nrow = 4)
    if (sum(choo3) != 0) {
      df_3way = vapply(which(choo3), FUN = function(combo) {
        df[, ways3[combo,1]]*df[, ways3[combo, 2]]*df[, ways3[combo, 3]]
      }, FUN.VALUE = rep(1,n))
      df_final = cbind(df_final, df_3way)
    } 
    
    
    if (way4) {
      df_4way = df[, 1]*df[, 2]*df[, 3]*df[, 4]
      df_final = cbind(df_final, df_4way)
    }
    
    if (truth) {
      N = 1e6
      df1 = vapply(1:d, FUN = function(x) funclist[[x]](Wbig[,x]), FUN.VALUE = rep(1, N))
      choo2 = x$choo2
      choo3 = x$choo3
      way4 = x$way4
      MT = x$MT
      funclist = x$funclist
      bin_coef = x$bin_coef
      # make number of interactions and which ones
      ways2 = matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), byrow = TRUE, nrow = 6)
      
      if (sum(choo2) != 0) {
        df_2way = vapply(which(choo2), FUN = function(combo) {
          df1[, ways2[combo,1]]*df1[, ways2[combo, 2]]
        }, FUN.VALUE = rep(1,N))
        df1_final = cbind(df1[,MT], df_2way)
      }
      
      # make number of 3 ways interactions
      ways3 = matrix(c(1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4), byrow = TRUE, nrow = 4)
      if (sum(choo3) != 0) {
        df_3way = vapply(which(choo3), FUN = function(combo) {
          df1[, ways3[combo,1]]*df1[, ways3[combo, 2]]*df1[, ways3[combo, 3]]
        }, FUN.VALUE = rep(1,N))
        df1_final = cbind(df1_final, df_3way)
      }
      
      
      if (way4) {
        df_4way = df1[, 1]*df1[, 2]*df1[, 3]*df1[, 4]
        df1_final = cbind(df1_final, df_4way)
      }
    }
    
    if (all(c(all(!choo2), all(!choo3), !way4))) {
      df_final = df
      df_final1 = df1_true
    }
    if (truth) return(list(df_final, df1_final)) else return(list(df_final, list()))
  })
  
  
  a = .7 
  coef_G = runif(ncol(dfs[[1]][[1]]), -a, a)
  PG_n  = gentmle2::truncate(plogis(dfs[[1]][[1]] %*% coef_G), .05)
  A = rbinom(n, 1, PG_n)
  mean(A)
  
  inters = rbinom(ncol(dfs[[2]][[1]]), 1, .5)
  dfinter_n = vapply(which(inters == 1), FUN = function(col) dfs[[2]][[1]][,col]*A, FUN.VALUE = rep(1, n))
  dfQn = cbind(dfs[[2]][[1]], dfinter_n)
  df1n = cbind(dfs[[2]][[1]], dfs[[2]][[1]][, inters])
  
  a = .3
  coef_Q = runif(ncol(dfQn), -a, a)

  # mean(Y)
  
  if (truth) {
    
    PG_true  = gentmle2::truncate(plogis(dfs[[1]][[2]] %*% coef_G), .05)
    Abig = rbinom(N, 1, PG_true)
    length(Abig)
    mean(Abig)
    
    
    dfinter_true = vapply(which(inters == 1), FUN = function(col) {
      dfs[[2]][[2]][,col]*Abig
      }, FUN.VALUE = rep(1, N))
    df1_true = cbind(dfs[[2]][[2]], dfs[[2]][[2]][, inters])
    dfQ_true = cbind(dfs[[2]][[2]], dfinter_true)
    
    BV0 = 0
    C = 1
    mm = 5
    jj = 1
    while (BV0 <= .025 & jj < 6) {
      coef_Q[(ncol(dfs[[2]][[1]]) + 1):ncol(dfQn)] = C*coef_Q[(ncol(dfs[[2]][[1]]) + 1):ncol(dfQn)]
      PQ_true1  = gentmle2::truncate(plogis(df1_true %*% coef_Q), .05)
      PQ_true0 = gentmle2::truncate(plogis(dfs[[2]][[2]] %*% coef_Q[1:ncol(dfs[[2]][[2]])]), .05)
      blip_true = PQ_true1 - PQ_true0
      ATE0 = mean(blip_true)
      BV0 = var(blip_true)
      C = C + .5
      jj = jj + 1
    }
  } 
  
  coef_Q[(ncol(dfs[[2]][[1]]) + 1):ncol(dfQn)] = coef_Q[(ncol(dfs[[2]][[1]]) + 1):ncol(dfQn)]
  PQ_n  = gentmle2::truncate(plogis(dfQn %*% coef_Q), .05)
  # hist(PQ_true0, breaks = 200)
  # BV0
  Y = rbinom(n, 1, PQ_n)
  
  if (!truth) {
    PQ_n1  = gentmle2::truncate(plogis(df1n %*% coef_Q), .05)
    PQ_n0 = gentmle2::truncate(plogis(dfs[[2]][[1]] %*% coef_Q[1:ncol(dfs[[2]][[1]])]), .05)
    blip_n = PQ_n1 - PQ_n0
    ATEn = mean(blip_n)
    BVn = var(blip_n)
  }
  if (truth) {
    return(list(BV0 = BV0, ATE0 = ATE0, W = W, A = A, Y = Y))
  } else {
    return(list(BVn = BVn, ATEn = ATEn, W = W, A = A, Y = Y))
  }
}

# A = get.info(1000, 4, FALSE) 
# 
# A$BVn
# A$ATEn
# mean(A$A)
# mean(A$Y)
# head(A$W)
# 
# A = get.info(1000, 4, TRUE) 
# A$BV0
# A$ATE0
# mean(A$A)
# mean(A$Y)
# head(A$W)

