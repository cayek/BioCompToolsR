#' for SFdS paper
#'
#'
#' @export
SFdS.runOnSimu <- function(K, d, rep, L, n0, n, L0, save.in.file, sigma = NULL, lambda = 1.0, max.iteration = 50) {

  # init
  seed = sample.int(.Machine$integer.max,1)
  # seed = runif(1)

  runOnSimu.res.rmse.n = data.table::data.table()
  runOnSimu.res.rmse.L = data.table::data.table()
  runOnSimu.res.normilized.residual.error = data.table::data.table()

  # on n
  for(i in n) {
    r = 1
    while(r <= rep) {
      # sink()
      cat("==n =",i,"rep =",r,"\n")

      #seed = seed + 1
      seed = sample.int(.Machine$integer.max,1)
      cat("seed :",seed,"\n")


      sink("/dev/null")
      set.seed(seed)
      data = sampleTESS2.3(i,L0,d,K,sample.dist.from.center.Q(0.1))
      W = heat.kernel.weight(data$coord,sigma)
      Lapl = graph.laplacian(W)

      set.seed(seed)
      projected.res = solverTESS3.projected(X = data$X,
                                            K = data$K,
                                            d = data$d,
                                            Lapl = Lapl,
                                            lambda = lambda,
                                            max.iteration = max.iteration)

      set.seed(seed)
      QP.res = solverTESS3.QP(X = data$X,
                                         K = data$K,
                                         d = data$d,
                                         Lapl = Lapl,
                                         lambda = lambda,
                                         max.iteration = max.iteration)
      sink()

      if (!QP.res$err) {
        runOnSimu.res.rmse.n = rbind(
          runOnSimu.res.rmse.n,
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$Q, data$Q),
                                  rep = r,
                                  n = K * i,
                                  L = L0,
                                  factor = "rmse(Q.QR,Q.true)"),
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$G, data$G),
                                  rep = r,
                                  n = K * i,
                                  L = L0,
                                  factor = "rmse(G.QR,G.true)"),
          data.table::data.table( rmse = rmse_withBestPermutation(projected.res$Q, data$Q),
                                  rep = r,
                                  n = K * i,
                                  L = L0,
                                  factor = "rmse(Q.proj,Q.true)"),
          data.table::data.table( rmse = rmse_withBestPermutation(projected.res$G, data$G),
                                  rep = r,
                                  n = K * i,
                                  L = L0,
                                  factor = "rmse(G.proj,G.true)")
        )

        runOnSimu.res.normilized.residual.error = rbind(
          runOnSimu.res.normilized.residual.error,
          data.table::data.table( normilized.residual.error =
                                    QP.res$normilized.residual.error,
                                  time = QP.res$time,
                                  rep = r,
                                  L = L0,
                                  n = K * i,
                                  algo = "QP",
                                  simu = paste("QP : n=", K*i, "|L=", L0,sep="")),
          data.table::data.table( normilized.residual.error =
                                    projected.res$normilized.residual.error,
                                  time = projected.res$time,
                                  rep = r,
                                  L = L0,
                                  n = K * i,
                                  algo = "proj",
                                  simu = paste("projected : n=", K*i, "|L=", L0,sep=""))
        )
        r = r + 1
      }
    }

  }

  # on L
  for(i in L) {
    r = 1
    while(r <= rep) {
      # sink()
      cat("==L =",i,"rep =",r,"\n")

      #seed = seed + 1
      seed = sample.int(.Machine$integer.max,1)
      cat("seed :",seed,"\n")

      sink("/dev/null")
      set.seed(seed)
      data = sampleTESS2.3(n0,i,d,K,sample.dist.from.center.Q(0.1))
      W = heat.kernel.weight(data$coord,sigma)
      Lapl = graph.laplacian(W)


      set.seed(seed)
      projected.res = solverTESS3.projected(X = data$X,
                                            K = data$K,
                                            d = data$d,
                                            Lapl = Lapl,
                                            lambda = lambda,
                                            max.iteration = max.iteration)

      set.seed(seed)
      QP.res = tryCatch({ solverTESS3.QP(X = data$X,
                                         K = data$K,
                                         d = data$d,
                                         Lapl = Lapl,
                                         lambda = lambda,
                                         max.iteration = max.iteration) },
                        error = function(e) {
                          sink()
                          cat("!! 1 error catched !! \n")
                          sink("/dev/null")
                        }, finally = {
                          err = TRUE
                        })
      sink()

      if (!QP.res$err) {
        runOnSimu.res.rmse.L = rbind(
          runOnSimu.res.rmse.L,
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$Q, data$Q),
                                  rep = r,
                                  n = K * n0,
                                  L = i,
                                  factor = "rmse(Q.QR,Q.true)"),
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$G, data$G),
                                  rep = r,
                                  n = K * n0,
                                  L = i,
                                  factor = "rmse(G.QR,G.true)"),
          data.table::data.table( rmse = rmse_withBestPermutation(projected.res$Q, data$Q),
                                  rep = r,
                                  n = K * n0,
                                  L = i,
                                  factor = "rmse(Q.proj,Q.true)"),
          data.table::data.table( rmse = rmse_withBestPermutation(projected.res$G, data$G),
                                  rep = r,
                                  n = K * n0,
                                  L = i,
                                  factor = "rmse(G.proj,G.true)")
        )

        runOnSimu.res.normilized.residual.error = rbind(
          runOnSimu.res.normilized.residual.error,
          data.table::data.table( normilized.residual.error =
                                    QP.res$normilized.residual.error,
                                  time = QP.res$time,
                                  rep = r,
                                  L = i,
                                  n = K * n0,
                                  algo = "QP",
                                  simu = paste("QP : L=", i, "|n=", K * n0,sep="")),
          data.table::data.table( normilized.residual.error =
                                    projected.res$normilized.residual.error,
                                  time = projected.res$time,
                                  rep = r,
                                  L = i,
                                  n = K * n0,
                                  algo = "proj",
                                  simu = paste("projected : L=", i, "|n=", K * n0,sep=""))
        )
        r = r + 1
      }

    }

  }

  # save
  runOnSimu.res = list(rmse.n=runOnSimu.res.rmse.n,
                       rmse.L=runOnSimu.res.rmse.L,
                       normilized.residual.error = runOnSimu.res.normilized.residual.error)
  save(runOnSimu.res, file = save.in.file)
  return(runOnSimu.res)
}

#' for SFdS paper
#'
#'
#' @export
SFdS.runOnTrueData <- function(K, rep, file.X, file.coord, save.in.file, sigma = NULL, lambda = 1.0, max.iteration = 50) {

  # init
  seed = sample.int(.Machine$integer.max,1)

  # read data
  X = LEA::read.geno(file.X)
  coord = tess3r::read.coord(file.coord)
  d = max(X)
  n = nrow(X)
  L = ncol(X)

  # compute graph
  W = heat.kernel.weight(coord,sigma)
  Lapl = graph.laplacian(W)

  # res data table
  runOnTrueData.res.rmse = data.table::data.table()
  runOnTrueData.res.normilized.residual.error = data.table::data.table()

  for(k in K) {
    r = 1
    while(r <= rep) {
      seed = sample.int(.Machine$integer.max,1)

      # sink()
      cat("==K =",k,"rep =",r,"\n")

      # sink("/dev/null")
      set.seed(seed)
      projected.res = solverTESS3.projected(X = X,
                                            K = k,
                                            d = d,
                                            Lapl = Lapl,
                                            lambda = lambda,
                                            max.iteration = max.iteration)

      set.seed(seed)
      QP.res = solverTESS3.QP(X = X,
                              K = k,
                              d = d,
                              Lapl = Lapl,
                              lambda = lambda,
                              max.iteration = max.iteration)
      # sink()
      if (!QP.res$err) {
        runOnTrueData.res.rmse = rbind(
          runOnTrueData.res.rmse,
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$Q, projected.res$Q),
                                  rep = r,
                                  K = k,
                                  n = n,
                                  L = L,
                                  factor = "rmse(Q.QR,Q.proj)"),
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$G, projected.res$G),
                                  rep = r,
                                  K = k,
                                  n = n,
                                  L = L,
                                  factor = "rmse(G.QR,G.proj)")
        )

        runOnTrueData.res.normilized.residual.error = rbind(
          runOnTrueData.res.normilized.residual.error,
          data.table::data.table( normilized.residual.error =
                                    QP.res$normilized.residual.error,
                                  rep = r,
                                  algo = paste("QP",sep=""),
                                  K = k,
                                  it = 1:max.iteration),
          data.table::data.table( normilized.residual.error =
                                    projected.res$normilized.residual.error,
                                  rep = r,
                                  algo = paste("projected",sep=""),
                                  K = k,
                                  it = 1:max.iteration)
        )
        r = r + 1
      }
    }
  }

  # save
  runOnTrueData.res = list(rmse=runOnTrueData.res.rmse,
                           normilized.residual.error = runOnTrueData.res.normilized.residual.error)
  save(runOnTrueData.res, file = save.in.file)
  return(runOnTrueData.res)

}

#' for SFdS paper
#'
#'
#' @export
SFdS.runOnTrueData.time <- function(K, file.X, file.coord, save.in.file, sigma = NULL, lambda = 1.0, max.iteration = 50) {

  rep = 1

  # times
  times = data.table::data.table()

  # init
  seed = sample.int(.Machine$integer.max,1)

  # read data
  X = LEA::read.geno(file.X)
  coord = tess3r::read.coord(file.coord)
  d = max(X)
  n = nrow(X)
  L = ncol(X)

  # compute graph
  W = heat.kernel.weight(coord,sigma)
  Lapl = graph.laplacian(W)

  # res data table
  runOnTrueData.res.rmse = data.table::data.table()
  runOnTrueData.res.normilized.residual.error = data.table::data.table()

  for(k in K) {
    r = 1
    while(r <= rep) {
      seed = sample.int(.Machine$integer.max,1)

      # sink()
      cat("==K =",k,"rep =",r,"\n")

      # sink("/dev/null")
      set.seed(seed)
      ptm <- proc.time()
      projected.res = solverTESS3.projected(X = X,
                                            K = k,
                                            d = d,
                                            Lapl = Lapl,
                                            lambda = lambda,
                                            max.iteration = max.iteration)
      t=proc.time() - ptm
      proj.time = t[1] + t[2] + t[4] + t[5]

      set.seed(seed)
      ptm <- proc.time()
      QP.res = solverTESS3.QP(X = X,
                              K = k,
                              d = d,
                              Lapl = Lapl,
                              lambda = lambda,
                              max.iteration = max.iteration)
      t=proc.time() - ptm
      QP.time = t[1] + t[2] + t[4] + t[5]

      # sink()
      if (!QP.res$err) {
        runOnTrueData.res.rmse = rbind(
          runOnTrueData.res.rmse,
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$Q, projected.res$Q),
                                  rep = r,
                                  K = k,
                                  n = n,
                                  L = L,
                                  factor = "rmse(Q.QR,Q.proj)"),
          data.table::data.table( rmse = rmse_withBestPermutation(QP.res$G, projected.res$G),
                                  rep = r,
                                  K = k,
                                  n = n,
                                  L = L,
                                  factor = "rmse(G.QR,G.proj)")
        )

        runOnTrueData.res.normilized.residual.error = rbind(
          runOnTrueData.res.normilized.residual.error,
          data.table::data.table( normilized.residual.error =
                                    QP.res$normilized.residual.error,
                                  rep = r,
                                  algo = paste("QP",sep=""),
                                  K = k,
                                  it = 1:max.iteration),
          data.table::data.table( normilized.residual.error =
                                    projected.res$normilized.residual.error,
                                  rep = r,
                                  algo = paste("projected",sep=""),
                                  K = k,
                                  it = 1:max.iteration)
        )

        times = rbind(times, data.table::data.table(K = k, time = proj.time / max.iteration, algo = "proj") )
        times = rbind(times, data.table::data.table(K = k, time = QP.time / max.iteration, algo = "QR") )

        r = r + 1
      }
    }
  }

  # save
  runOnTrueData.res = list(times = times, rmse=runOnTrueData.res.rmse,
                           normilized.residual.error = runOnTrueData.res.normilized.residual.error)
  save(runOnTrueData.res, file = save.in.file)
  return(runOnTrueData.res)

}

