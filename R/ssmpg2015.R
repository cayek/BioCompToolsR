#' Compute true power
#'
#'
cleanAnswer <- function( dirtyCsv ) {
  aux=dirtyCsv
  aux$Answer = lapply(strsplit(as.character(dirtyCsv$Answer),","), as.numeric)
  return(aux)
}


#' compute precision and recall
#'
#'
#' @export
fdrpreCalc = function(csvfile,ground_truth){
  aux = csvfile
  aux$Answer = lapply(strsplit(as.character(csvfile$Answer),","),as.numeric)
  N <- length(aux$Answer)
  FDR <- NULL
  Power <- NULL
  for (n in 1:N){
    answer <- unlist(aux$Answer[n])
    if (length(answer)==0){
      FDR <- c(FDR,0)
      Power <- c(Power,0)
    } else {
      FDR <- c(FDR,sum(!(answer %in% ground_truth))/length(answer))
      Power <- c(Power,sum((answer %in% ground_truth))/length(ground_truth))
    }
  }
  return(data.frame(csvfile,FDR,Power))
}



#' Download contest result
#'
#'
#' @export
ddl.result <- function( case.nb, ground_truth ) {

  url = paste("https://ssmpg-challenge.imag.fr/challenge/",
              case.nb,"/api/csv_generate",sep="" )

  tmpFile <- tempfile()
  download.file(url, destfile = tmpFile, method = "curl")
  result <- read.csv(tmpFile, sep = "\t",
                     colClasses=c('character',
                                  'character',
                                  'numeric',
                                  'character',
                                  'factor',
                                  'factor',
                                  'character',
                                  'factor',
                                  'character'))

  result_withfdrtp = fdrpreCalc(result,ground_truth)
  result_withfdrtp$Date = as.POSIXct(result_withfdrtp$Date, "%Y-%m-%d %H:%M:%S+00:00")
  return(result_withfdrtp)
}

#' keep only best submision by method and the
#'
#'
#' @export
filterbest <- function( result ) {

  result.best = result
  indices = 1:nrow(result)
  toremove = c()

  best_methodxuserxenv = merge( aggregate(Date ~Method + User + EnvVar ,
                                data = result,
                                max ), result)
  best_userxenv = merge( aggregate(Score ~ User + EnvVar ,
                                          data = best_methodxuserxenv,
                                          max ), best_methodxuserxenv)

  return(best_userxenv)
}



#' give the winner
#'
#'
#' @export
winner <- function( result ) {

  winnerresult = sort(result)
  winnerresult = winnerresult[winnerresult$Score == max(winnerresult$Score),]
  return(winnerresult)
}

#' sort
#'
#'
#' @export
sortWinner <- function( result,
                  score = function(fdr,precision) {
                    return(2*((1-fdr)*precision)/((1-fdr)+precision))}, rankFun = sum ) {

  winnerresult = filterbest(result)
  # remove method if there is several max
  winnerresult = aggregate(Score ~ User+ EnvVar ,
                           data = winnerresult,
                           mean )
  winnerresult = aggregate(Score ~ User ,
                           data = winnerresult,
                           rankFun )
  return(winnerresult[order(winnerresult$Score),])
}

#' cheater
#'
#'
#' @export
duplitedRowDetection <- function( result ) {

  return(result[duplicated(result[,c("User","Method","Score","Answer","Desc")]),])

}


#' noise
#'
#'
#' @export
removeNoise <- function( result ) {
  # remove duplicate
  result = result[!duplicated(result[,c("User","Method","Score","Answer","Desc")]),]
  return(result)

}

