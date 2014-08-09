plotVar.cia <-
function(x, var=NA, axes=1:2,
         var.col="red", # the length either 1 or length(var)
         var.lab=FALSE, # T or F
         bg.var.col="gray", # the length either 1 or 2
         nlab=0,
         sepID.data=NULL, # NULL, 1, 2, 1:2
         sepID.sep="_", ...) 
  {
  
  coin <- x
  if (length(bg.var.col) == 1)
    bg.var.col <- rep(bg.var.col, 2)
  if (length(var.col) == 1)
    var.col <- rep(var.col, length(var))
  
  df1 <- coin$coinertia$co[axes]
  df2 <- coin$coinertia$li[axes]
  
  ndf1 <- rownames(df1)
  if (1 %in% sepID.data) 
    ndf1 <- sapply(lapply(strsplit(ndf1, split=sepID.sep), 
                          function(x) x[-length(x)]), paste, collapse=sepID.sep)
  ndf2 <- rownames(df2)
  if (2 %in% sepID.data) 
    ndf2 <- sapply(lapply(strsplit(ndf2, split=sepID.sep), 
                          function(x) x[-length(x)]), paste, collapse=sepID.sep)
  
  layout(1:2)
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  
  plotgenes(df1, axis1=axes[1], axis2=axes[2], nlab=nlab, genelabels=ndf1, colpoints=bg.var.col[1])
  ind <- ndf1 %in% var
  if (any(ind)) {
    points(df1[ind, ], pch=20, col=var.col[na.omit(match(ndf1, var))])
    if (var.lab)
      text(df1[ind, 1], df1[ind, 2], ndf1[ind])
  }
  legend(x="bottomleft", bty="n", legend="data1", x.intersp=-.5)
  
  plotgenes(df1, axis1=axes[1], axis2=axes[2], nlab=nlab, genelabels=ndf1, colpoints=bg.var.col[1])
  ind <- ndf2 %in% var
  if (any(ind)) {
    points(df2[ind, ], pch=20, col=var.col[na.omit(match(ndf2, var))])
    if (var.lab)
      text(df2[ind, 1], df2[ind, 2], ndf2[ind])
  }
  legend(x="bottomleft", bty="n", legend="data2", x.intersp=-.5)
  
  rep <- cbind(var, var %in% ndf1, var %in% ndf2)
  colnames(rep) <- c("variables", "data1", "data2")
  rep <- as.data.frame(rep)
  if (!is.na(var[1]))
    return(rep)
}
