selectVar.mcia <-
function(x, axis1=1, axis2=2,
         a1.lim=c(-Inf, Inf), a2.lim=c(-Inf, Inf), 
         sepID.data=NULL, sepID.sep="_", ...) # char or numeric
  {
  mcoin <- x
  if (!inherits(mcoin, "mcia"))
    stop("mcia object expected, please run mcia first")
  
  lab <- mcoin$mcoa$TC$"T"
  coord <- mcoin$mcoa$Tco
  coa <- names(mcoin$coa)
  
  if (is.list(a1.lim)) {
    if (length(a1.lim) != length(coa)) 
      stop ("length of a1.lim doesn't equals the number of datasets")
    sapply(a1.lim, .checklim)
  } else if (is.numeric(a1.lim)) {
    .checklim(a1.lim)    
    temp <- lapply(coa, function(x) a1.lim)
    a1.lim <- temp
  } else 
    stop("either numeric of a list of numeric required for a1.lim")
  
  if (is.list(a2.lim)) {
    if (length(a2.lim) != length(coa)) 
      stop ("length of a2.lim doesn't equals the number of datasets")
    sapply(a2.lim, .checklim)
  } else if (is.numeric(a2.lim)) {
    .checklim(a2.lim)    
    temp <- lapply(coa, function(x) a2.lim)
    a2.lim <- temp
  } else 
    stop("either numeric of a list of numeric required for a1.lim")
  
  j <- FALSE
  rowsel <- list()
  for (i in 1 : length(coa)) {
    sel <- lab %in% levels(lab)[i] 
    a1bottom <- coord[, axis1] > a1.lim[[i]][1]
    a1up <- coord[, axis1] < a1.lim[[i]][2]
    a1 <- a1bottom & a1up & sel
    
    a2bottom <- coord[, axis2] > a2.lim[[i]][1]
    a2up <- coord[, axis2] < a2.lim[[i]][2]
    a2 <- a2bottom & a2up & sel
    
    var <- rownames(coord) [a1 & a2]
    if (length(var) > 0)
      var <- sapply(strsplit(var, "\\."), function(x) x[1])    
    rowsel[[i]] <- var
    j <- j | any(a1 & a2)
  }
  
  if (!is.null(sepID.data))
    for (i in sepID.data) {
      if (length(rowsel[[i]]) > 0)
        rowsel[[i]] <- unique(sapply(lapply(strsplit(rowsel[[i]], split=sepID.sep), function(x) x[-length(x)]), paste, collapse=sepID.sep))

    }
    
  
  if (!j) {
    mat <- NULL
    print("NO variables selected in this region")
  } else {
    var <- unique(unlist(rowsel))
    mat <- sapply(rowsel, function(x) var %in% x)
    if (is.vector(mat))
      mat <- matrix(mat, 1)
    colnames(mat) <- coa
    # rownames(mat) <- allvar
    stat <- apply(mat, 1, function(x) length(which(x)))
    mat <- cbind(var, data.frame(mat), stat)
    mat$var <- as.character(mat$var)
  }
  return(mat)
}
