plotVar.mcia <-
function(x, var=NA, axes=1:2, 
         var.col="red", # the length either 1 or length(var)
         var.lab=FALSE, # T or F
         bg.var.col="gray", # the length either 1 or length(df)
         nlab=0, 
         sepID.data=NULL, 
         sepID.sep="_", 
         df=NA, # either name of data.frame or numeric
         layout=NA, ...) 
  
  {
  mcoin <- x
  rows <- mcoin$mcoa$Tco[axes]
  row.lab <- mcoin$mcoa$TC$"T"
  datasets <- names(mcoin$coa)
  
  if (is.na(df[1]))
    df <- 1:length(datasets) else
      if (!all(df %in% 1:length(datasets)) & !all(df %in% datasets))
        stop("undefined data.frame selected")
  
  if (!inherits(mcoin, "mcia"))
    stop("mcia object expected for mcoin")
  if (!is.null(sepID.data[1]) && any(grep(sepID.sep, var)))
    warning("be careful: sepID.separator present in variable names")
  if (length(axes) != 2)
    stop("you have to select (only) 2 axis")
  if (length(var.col) != 1 & length(var.col) != length(var))
    stop("the length of var.col could only be either 1 or length(var)")
  if (length(bg.var.col) != 1 & length(bg.var.col) != length(df))
    stop("the length of bg.var.col could only be either 1 or length(df)")
  
  if (length(bg.var.col) == 1)
    bg.var.col <- rep(bg.var.col, length(df))
  if (length(var.col) == 1)
    var.col <- rep(var.col, length(var))
  
  df.list <- lapply(levels(row.lab), 
                    function(x) return(rows[row.lab %in% x, ]))
  
  
  n <- length(df)
  n <- ceiling(n/2)*2
  
  rep <- var
  if (is.matrix(layout))
    layout(layout) else
      if (is.na(layout)) {
        if (length(df) == 1)
          layout(1) else
            if (length(df) == 2)
              layout(t(t(1:2))) else
                if (length(df) == 3)
                  layout(t(t(1:3))) else
                    if (length(df) > 3)
                      layout(matrix(1:n, n/2, byrow=T))
      }
  
  par(mar=c(0.1, .1, .1, .1))
  for (i in df) {
    where <- match(i, df)
    idf <- df.list[[i]]
    ns <- rownames(idf)
    ns <- sapply(lapply(strsplit(ns, split="\\."), 
                        function(x) x[-length(x)]), paste, collapse=".")
    plotgenes(idf, axis1=1, axis2=2, nlab=nlab, genelabels=ns, colpoints=bg.var.col)
    
    if (!is.null(sepID.data))
      if (i %in% sepID.data)
        ns <- sapply(lapply(strsplit(ns, split=sepID.sep), 
                            function(x) x[-length(x)]), paste, collapse=sepID.sep)
    ind <- ns %in% var
    if (any(ind)) {
      points(idf[ind, ], pch=20, col=var.col[na.omit(match(ns, var))])
      if (var.lab)
        text(idf[ind, 1], idf[ind, 2], ns[ind])
    }
    legend(x="bottomleft", bty="n", legend=datasets[i], x.intersp=-.5)
    rep <- cbind(rep, var %in% ns)
  }
  colnames(rep) <- c("Variables", datasets[df])
  rep <- as.data.frame(rep)
  if (!is.na(var[1]))
    return(rep)
}
