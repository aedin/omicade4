plot.mcia <-
function(x, axes=1:2, sample.lab=TRUE, sample.legend=FALSE, sample.color=1, phenovec=NULL, df.color=1, df.pch=NA, gene.nlab=0, ...) 
  {
  mcoin <- x
  if (!inherits(mcoin, "mcia"))
    stop("mcia object expected, please run mcia first")
  
  ndata <- length(mcoin$coa)
  eig <- mcoin$mcoa$pseudoeig
  cov2 <- mcoin$mcoa$cov2[, axes] #pseueig
  

  
  if (!length(sample.color) %in% c(1, nrow(mcoin$mcoa$SynVar)))
    stop("length of sample.color should be either 1 or # of samples")
  if (!length(df.color) %in% c(1, ndata))
    stop("length of df.color should be either 1 or # of samples")
  
  if (is.na(df.pch[1])) {
    pch <- c(16, 17, 15, 18, 1, 2, 0, 5)
    if (ndata > 8) {
      pch <- rep(pch, ceiling(ndata/8))[1:ndata]
      warning("more than 8 datasets in mcia, df.pch is recycled used")
    } else {
      pch <- pch[1:ndata]
    }
  } else {
    if (length(df.pch) != ndata)
      stop("the length of df.pch should be the same with datasets in mcia, recycled use df.pch is not allowed")
    pch <- df.pch
  }
         
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  splot.sample.mcia(mcoin, axis1=axes[1], axis2=axes[2], 
                    col=sample.color, pch=pch, 
                    sample.lab=sample.lab, 
                    legend=sample.legend, 
                    sub="sample space",
                    phenovec=phenovec)
  
  splot.mol.mcia(mcoin, axis1=axes[1], axis2=axes[2], col=df.color, pch=pch, gene.nlab=gene.nlab)
  
  #plot c(2, 1)
#   par(mar=c(3, 3, 3, 3))
#   barplot(eig, col="black",names.arg=paste("eig", 1:length(eig)))
#   legend(x="topright", pch=pch, col=df.color, legend=names(mcoin$coa), box.col="white")
  
  nkeig <- ncol(mcoin$mcoa$Tli)
  eig <- mcoin$mcoa$pseudoeig
  proe <- eig/sum(eig)
  
  par(mar=c(3, 4, 1, 4))
  neig <- length(eig)
  po <- barplot(eig, plot=FALSE)
  barplot(eig,names.arg= 1:neig, col=c(rep("cyan", nkeig), rep("gray", neig-nkeig)), xlab="", ylab="Eigen Value")
  
  par(new=T) 
  plot(cbind(po, proe), frame.plot=FALSE, pch=20, axes=FALSE, xlab="Eigen Vector", ylab="")
  points(x=po, y=proe, pch=20)
  lines(x=po, y=proe)
  axis(side=4)
  mtext(side=4, text="Percentage", line=2.5, cex=0.8)
  legend(x="topright", pch=pch, col=df.color, legend=names(mcoin$coa), box.col="white")
  
  #plot c(2, 2)
  par(mar=c(4.5, 4.5, 0.5, 0.5))
  plot(cov2, pch=".", main="", axes=TRUE, col=NA, 
       xlab=paste("pseudoeig", axes[1]), ylab=paste("pseudoeig", axes[2]))
  scatterutil.grid(0)
  points(cov2, pch=pch, col=df.color)
}
