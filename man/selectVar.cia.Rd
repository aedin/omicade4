\name{selectVar.cia}
\alias{selectVar.cia}
\title{
  Selecting variables from result of CIA
}
\description{
  To select variables in CIA
  variable space, the function is called by \code{selectVar}.
}
\usage{
\method{selectVar}{cia}(x, axis1 = 1, axis2 = 2, 
  df1.a1.lim = c(-Inf, Inf), df1.a2.lim = c(-Inf, Inf), 
  df2.a1.lim = df1.a1.lim, df2.a2.lim = df1.a2.lim, 
  sepID.data = NULL, sepID.sep = "_", ...)
}

\arguments{
  \item{x}{
    The result returned by \code{\link[made4]{cia}}
}
  \item{axis1}{
    Integer, the column number for the x-axis. The default is 1.
}
  \item{axis2}{
    Integer, the column number for the y-axis. The default is 2.
}
  \item{df1.a1.lim}{
    A vector containing 2 numbers indicating the range of \bold{X} axis 
    of selected on the \bold{1st} data.frame. The first value limiting the 
    lower boundary, the second value limiting the upper boundary.
}
  \item{df1.a2.lim}{
    The range of \bold{Y} axis of selected on the \bold{1st} datasets. 
}
  \item{df2.a1.lim}{
    The range of \bold{X} axis of selected on the \bold{2nd} dataset. 
}
  \item{df2.a2.lim}{
    The range of \bold{Y} axis of selected on the \bold{2nd} dataset. 
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
  \code{cia} projecting variables from different datasets
  to a two dimensional space. This function supplies a method selecting 
  variables according to the co-ordinates of variables
}
\value{
  Returns a data.frame describing which variables are presented
  on which data.frame within the limited region(s).
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
  See Also as \code{\link{selectVar.mcia}}
}

