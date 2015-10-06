#' Plot pc1 and pc2
#'
#' @export
plotPc <- function( genotype, scale = FALSE, xpcindex = 1, ypcindex = 2, ... ) {

  pc = prcomp(genotype, scale = scale) # Remark when setting T : there are some outliers but if we set F they disappear !!
  ggplot2::qplot( x = pc$x[,xpcindex], y = pc$x[,ypcindex], ...) +
    ggplot2::theme(legend.position="none")

}


#' Plot folded allele frequency spectrum
#'
#' @param genotype A matrix of genotype of size nbIndiv x nbLocus
#'
#' @export
plotSNPfoldedspectrum <- function( genotype ) {
  spectrum = apply( data$genotype, MARGIN = 2, FUN = function(x) min(sum(x==1),sum(x==0)))
  plot(table(spectrum))
}


#' The allele frequency spectrum
#'
#' @param genotype A matrix of genotype of size nbIndiv x nbLocus
#' @param allel 0 or 1 for allele which you want the spectrum
#'
#' @export
plotSNPspectrum <- function( genotype, allele = 1 ) {
  spectrum = apply( data$genotype, MARGIN = 2, FUN = function(x) sum(x==allele))
  plot(table(spectrum))
}


#' Manhattan plot
#'
#' @export
manhattanPlot <- function( p.values ) {

  ggplot( data = data.frame(p.values = -log(p.values), index = 1:length(p.values) ), aes( x = index, y = p.values ) ) +
    geom_point(colour = rainbow(1))

}



#' Multiple graphs on one page
#'
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

    # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

