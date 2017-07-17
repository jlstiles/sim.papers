#' @export
coverage.check = function(data, truth, ind) {
  ans = vapply(ind,FUN = function(x){
    covs = data[,x+1]<=var0&data[,x+2]>=var0
    mean(covs)
  }, FUN.VALUE = 1)
  names(ans) = colnames(data)[ind]
  return(ans)
}

#' @export
coverage.checksimul = function(data, truth, ind) {
  covs = data[,ind[1]+1]<=var0&data[,ind[1]+2]>=var0
  covs1 = data[,ind[2]+1]<=var0&data[,ind[2]+2]>=var0
  mean(covs*covs1)
}
