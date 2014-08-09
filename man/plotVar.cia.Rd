\name{plotVar.cia}
\alias{plotVar.cia}
\title{
  Plot variable space of result from \code{\link[made4:cia]{cia}}
}
\description{
  Plot variable space of \code{\link[made4:cia]{mcia}} 
  and visualize selected variables across datasets.
}
\usage{
\method{plotVar}{cia}(x, var = NA, axes = 1:2, 
  var.col = "red", var.lab = FALSE, bg.var.col = "gray", 
  nlab = 0, sepID.data = NULL, sepID.sep = "_", ...)
}
\arguments{
  \item{x}{
    An object of class \code{cia}
}
  \item{var}{
    see \code{\link{plotVar}}
}
  \item{axes}{
    see \code{\link{plotVar}}
}
  \item{var.col}{
    see \code{\link{plotVar}}
}
  \item{var.lab}{
    see \code{\link{plotVar}}
}
  \item{bg.var.col}{
    see \code{\link{plotVar}}
}
  \item{nlab}{
    see \code{\link{plotVar}}
}
  \item{sepID.data}{
    see \code{\link{plotVar}}
}
  \item{sepID.sep}{
    see \code{\link{plotVar}}
}
  \item{\dots}{
    Other arguments
}
}

\value{
  If \code{var} is not NA, a data frame is return, with rows for variables of
  interest and columns of logical value indicating which \code{data.frame}s contains which
  variables.
}
%\references{
%  Meng C, Kuster B, Culhane AC, and Gholami AM. A multivariate 
%  approach to the integration of multi-omics datasets. 
%  BMC bioinformatics (under review)
%}
\author{
  Chen Meng
}

\seealso{
  See Also as \code{\link{plotVar.mcia}}
}

