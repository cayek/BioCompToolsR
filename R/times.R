#' compute run time
#'
#'
#' @export
computeRuntime <- function(fun,rep = 1) {

  res = 1:rep
  for(i in 1:rep) {
    ptm <- proc.time()
    fun()
    t=proc.time() - ptm
    res[i]=t[4] + t[5]+ t[1] + t[2]
  }
  return(res)
}
