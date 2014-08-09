\name{selectVar.mcia}
\alias{selectVar.mcia}
\title{
  Selecting variables from result of MCIA
}
\description{
  The selection of variables based on co-ordinates of MCIA
  variable space. The function is called by \code{selectVar}
}
\usage{
\method{selectVar}{mcia}(x, axis1 = 1, axis2 = 2, 
  a1.lim = c(-Inf, Inf), a2.lim = c(-Inf, Inf), 
  sepID.data = NULL, sepID.sep = "_", ...)
}

\arguments{
  \item{x}{
    An object of class \code{mcia}, the result returned by \code{\link{mcia}}.
}
  \item{axis1}{
    Integer, the column number for the x-axis. The default is 1.
}
  \item{axis2}{
    Integer, the column number for the y-axis. The default is 2.
}
  \item{a1.lim}{
    The limited range of x-axis of selected. It could be either a vector (containing 
    2 numbers, the first value limiting the lower boundary, the second 
    value limiting the upper boundary) or a list of vectors, each of which 
    contains two number. If it is a \code{list}, the length of the \code{list} should be 
    the same with number of \code{data.frame}s in \code{mcia}.
}
  \item{a2.lim}{
    The limited range of y-axis.
}
  \item{sepID.data}{
    See \code{\link{plotVar.mcia}}
}
  \item{sepID.sep}{
    See \code{\link{plotVar.mcia}}
}
  \item{\dots}{
    Other arguments
}
}
\details{
  \code{mcia} projecting variables (genes) from different datasets
  to a lower dimensional space. This function supplies a method selecting 
  variables according to the co-ordinates of variables.
}
\value{
  Returns a data.frame describing which variables are presented
  on which data.frames within the limited region(s).
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
  See Also as \code{\link{selectVar.cia}},  \code{\link{selectVar}}
}
\examples{
  data(NCI60_4arrays)
  mcoin <- mcia(NCI60_4arrays)
  selectVar(mcoin, a1.lim=c(1, Inf))
}

