#' compute the root mean square error
#'
#'
#' @export
rmse <- function( Q1, Q2) {
  return(sqrt(mean( (Q1 - Q2)^2)))
}

#' Compute the root mean square error between matrix by permuting
#' matrix column such that we have the best rmse.
#'
#'
#'
#' @export
rmse_withBestPermutation <- function( Q1, Q2) {

  library("permute")

  aux = rmse(Q1,Q2)

  K = dim(Q1)[2]
  #Because of R !!!!
  if( K == 2 ) {
    perms = matrix(c(2,1),nrow = 1,ncol=2)
  } else {
    perms = allPerms(K)
  }
  for(i in 1:(dim(perms)[1])) {

    aux1 = rmse(Q1,Q2[,perms[i,]])
    if(aux1 < aux) {
      aux = aux1
    }

  }
  return(aux)
}

#' Compute the root mean square error between matrix by permuting
#' matrix column such that we have the best rmse. Warning : Greedy algorithm !
#'
#'
#'
#' @export
rmse_withBestPermutation.greedy <- function( Q1, Q2) {

  # find perm
  perm = 1:ncol(Q1)

  taken = c()

  for(i in 1:ncol(Q1)) {
    min = .Machine$double.xmax
    if(i>1) {
      taken=c(taken,perm[i-1])
    }
    for(j in 1:ncol(Q1)) {
      aux = rmse(Q1[,i],Q2[,j])
      if(aux < min && !(j %in% taken)) {
        perm[i] = j
        min = aux
      }
    }
  }

  if(sum(duplicated(perm)) != 0){
    stop("this is not a permutation")
  }
  return(rmse(Q1,Q2[,perm]))
}
