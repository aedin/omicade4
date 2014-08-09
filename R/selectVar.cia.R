selectVar.cia <-
function(x, axis1=1, axis2=2, 
         df1.a1.lim=c(-Inf, Inf), df1.a2.lim=c(-Inf, Inf), 
         df2.a1.lim=df1.a1.lim, df2.a2.lim=df1.a2.lim, 
         sepID.data=NULL, sepID.sep="_", ...) # only numeric 
    {
    coin <- x
    if (!inherits(coin, "cia"))
      stop("cia object expected, please run cia first")
    
    .checklim(df1.a1.lim)
    .checklim(df1.a2.lim)
    .checklim(df2.a1.lim)
    .checklim(df2.a2.lim)
    
    
    
    df1 <- coin$coinertia$co
    df2 <- coin$coinertia$li
    
    df1.var <- rownames(df1)
    df2.var <- rownames(df2)
    
    if (length(df1.var) == 0 & length(df2.var) == 0) {
      print("NO variables selected in this region")
      mat <- NULL
    } else {
      if (! is.null(sepID.data)) {
        if (! is.numeric(sepID.data)) 
          stop("sepID.data should be 1, 2 or 1:2, which indicates the first, second or both of
               datasets are sepID.datasets repectively")
        if (sepID.data == 1) 
          df1.var <- sapply(lapply(strsplit(df1.var, split=sepID.sep), function(x) x[-length(x)]), paste, collapse=sepID.sep)
        if (sepID.data == 2) 
          df2.var <- sapply(lapply(strsplit(df2.var, split=sepID.sep), function(x) x[-length(x)]), paste, collapse=sepID.sep)
      }
      
      a1.1 <- df1[, axis1] > df1.a1.lim[1]
      a1.2 <- df1[, axis1] < df1.a1.lim[2]
      a1.3 <- df1[, axis2] > df1.a2.lim[1]
      a1.4 <- df1[, axis2] < df1.a2.lim[2]
      
      a2.1 <- df2[, axis1] > df2.a1.lim[1]
      a2.2 <- df2[, axis1] < df2.a1.lim[2]
      a2.3 <- df2[, axis2] > df2.a2.lim[1]
      a2.4 <- df2[, axis2] < df2.a2.lim[2]
      
      smdf1 <- df1.var[a1.1 & a1.2 & a1.3 & a1.4]
      smdf2 <- df2.var[a2.1 & a2.2 & a2.3 & a2.4]
      
      var <- unique(c(smdf1, smdf2))
      p1 <- var %in% smdf1
      p2 <- var %in% smdf2
      stat <- as.numeric(p1 & p2) + 1
      mat <- cbind(var, p1, p2, stat)
      mat <- as.data.frame(mat)
      colnames(mat) <- c("variables", "coa1", "coa2", "stat")
    }
    return(mat)
  }
