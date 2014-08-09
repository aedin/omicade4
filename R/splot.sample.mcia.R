splot.sample.mcia <-
function(x, axis1=1, axis2=2, 
                              col=1, pch=20, 
                              sample.lab=TRUE, 
                              legend=TRUE, 
                              sub="",
                              phenovec=NULL) {
  
  # plot matched samples from mcia
  if (!inherits(x, "mcia"))
    stop("x should be an object of class mcia!")
  if (!is.null(phenovec))
    phenovec <- as.factor(phenovec)
  ndata <- length(x$coa)
  dfxy <- x$mcoa$Tl1
  syn <- x$mcoa$SynVar
  pl <- pch
  cl <- col
  
  if (is.logical(legend)) {
      if (legend && is.null(phenovec) && length(col) == 1)
    warning("legend=TRUE is ignored since neither col nor phenovec is used to distinguish
      observations.")
  }
  if (!is.null(phenovec) && length(phenovec) != nrow(syn))
    stop("the length of phenovec should be the same with # of samples")
  if (!axis1 %in% 1:ncol(dfxy))
    stop("Uncharacterized axis selected by axis1")
  if (!axis2 %in% 1:ncol(dfxy))
    stop("Uncharacterized axis selected by axis1")
  if (!length(col) %in% c(1, nrow(syn), length(levels(phenovec))))
    stop("length of col should be either 1 or # of samples")
  if (!length(pch) %in% c(1, ndata))
    stop("length of pch should be either 1 or # of data frame")
  
  sync <- c()
  for (i in 1:(ndata)) {
    sync <- rbind(sync, syn)
  }
  c <- x$mcoa$TL$"T"
  
  if (length(col) == 1) {
    if (is.null(phenovec))
      col <- rep(col, length(c)) else
        col <- c(phenovec)
  } else if (length(col) == length(levels(phenovec))) {
    if (!is.null(phenovec))
      col <- col[c(phenovec)]
  } else
    col <- rep(col, ndata)
  
  if (length(pch) == 1)
    pch <- rep(pch, length(c)) else
      pch <- rep(pch, table(c))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  coo <- scatterutil.base(dfxy = dfxy, xax = axis1, yax = axis2, sub = sub, 
                          xlim = NULL, ylim = NULL, grid = TRUE, 
                          addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
                          origin = c(0, 0), csub = 1.25, possub = "bottomleft", 
                          pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE)
  
  points(dfxy[, c(axis1, axis2)], pch=pch, col=col)
  segments(sync[, axis1], sync[, axis2], dfxy[, axis1], dfxy[, axis2], col=col)
  
  if (sample.lab) {
    lab <- rownames(syn)
    text(x=syn[, axis1], y=syn[, axis2], lab)
  }
  
  if (is.logical(legend)) {
    if (legend && any(c(length(cl) != 1, !is.null(phenovec)))) {
      
      cle <- c()
      if (length(cl) == nrow(syn)) {
        cl <- cl
        cle <- rownames(syn)
      } else if (length(cl) == 1 && !is.null(phenovec)) {
        cl <- sort(unique(c(phenovec)))
        cle <- levels(phenovec)
      } else if (length(cl) == length(levels(phenovec))) {
        cl <- cl
        cle <- levels(phenovec)
      }
      pch.i <- 20
      col.i <- cl
      le.i <- cle
      
      legend("topleft", fill=FALSE, col=col.i, pch=pch.i, legend=le.i, border=F, bty="n")
    }
  } else if (is.character(legend) & length(legend) == 1) {
    eval(parse(text=legend))
  } else 
    stop("unknown legend setting")
  box()
}
