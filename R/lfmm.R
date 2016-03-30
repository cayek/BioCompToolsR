#' @export
lfmm.wrapper <- function(genome, env,...) {

  tmpDir = "~/PatatorHomeDir/tmp/"

  tmp = tempfile(tmpdir = tmpDir)
  temp.lfmm = paste(tmp,".lfmm",sep="")
  temp.env = paste(tmp,".env",sep="")
  LEA::write.env(env, temp.env)
  LEA::write.lfmm(genome, temp.lfmm)
  lfmm.obj = LEA::lfmm(temp.lfmm,temp.env,...)
  file.remove(temp.lfmm)
  file.remove(temp.env)
  return(lfmm.obj)

}

#' @export
onerun.lfmm.wrapper <- function(genome, env,parameter) {

  if(is.null(parameter$K)) {
    stop("Parameter K is mandatory")
  }

  tmpDir = "~/PatatorHomeDir/tmp/"

  tmp = tempfile(tmpdir = tmpDir)
  temp.lfmm = paste(tmp,".lfmm",sep="")
  temp.env = paste(tmp,".env",sep="")
  LEA::write.env(env, temp.env)
  LEA::write.lfmm(genome, temp.lfmm)

  lfmm.obj = LEA::lfmm(temp.lfmm,temp.env,K = parameter$K, project = "new")
  pvalue = as.vector(LEA::p.values(lfmm.obj, K = parameter$K, run = 1))
  zs = as.vector(LEA::z.scores(lfmm.obj, K=parameter$K, run = 1))

  # remove tmp file
  file.remove(temp.lfmm)
  file.remove(temp.env)
  return(list(pvalue = pvalue, zs = zs))

}


#' @export
combined.lfmm.wrapper <- function(genome, env,parameter) {

  if(is.null(parameter$K)) {
    stop("Parameter K is mandatory")
  }

  tmpDir = "~/PatatorHomeDir/tmp/"
  # create tmp file
  tmp = tempfile(tmpdir = tmpDir)
  temp.lfmm = paste(tmp,".lfmm",sep="")
  temp.env = paste(tmp,".env",sep="")
  LEA::write.env(env, temp.env)
  LEA::write.lfmm(genome, temp.lfmm)
  lfmm.obj = LEA::lfmm(temp.lfmm,temp.env,K=parameter$K,rep = 5, project = "new")

  #Record z-scores from the 5 runs in the zs matrix
  zs = LEA::z.scores(lfmm.obj, K=parameter$K)

  #Combine z-scores using the median
  zs.median = apply(zs, MARGIN = 1, median)
  #Compute the GIF
  lambda = median(zs.median^2)/0.456
  lambda
  # compute adjusted p-values from the combined z-scores
  adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)

  # remove file
  file.remove(temp.lfmm)
  file.remove(temp.env)
  return(list(pvalue = adj.p.values, zs.median = zs.median))

}
