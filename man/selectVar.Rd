\name{selectVar}
\alias{selectVar}
\title{
  Selecting variables (genes) from result of MCIA or CIA according to co-ordinates
}
\description{
  The user level function calls \code{selectVar.mcia} or 
  \code{selectVar.cia}. Function \code{cia} or \code{mcia} projects variables (genes) 
  from different datasets
  to a 2 dimensional space. This function supplies a method selecting 
  variables (genes) according to the coordinates of variables
}
\usage{
  selectVar(x, axis1 = 1, axis2 = 2, ...)
}
\arguments{
  \item{x}{
    An object of class \code{cia} or \code{mcia}, the result returned by 
    \code{\link[made4]{cia}} or \code{\link{mcia}} respectively.
}
  \item{axis1}{
    Integer, the column number for the x-axis. The default is 1.
}
  \item{axis2}{
    Integer, the column number for the y-axis. The default is 2.
}
  \item{\dots}{
    Other arguments
}
}

\value{
  Returns a data.frame describing which variables (genes) are presented
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
  See Also as \code{\link{selectVar.mcia}}, \code{\link{selectVar.cia}}
}


\examples{

data(NCI60_4arrays)
mcoin <- mcia(NCI60_4arrays)
selectVar(mcoin, a1.lim=c(2, Inf), a2.lim=c(-Inf, Inf))

# an example for the usage of sepID.data and sepID.sep
nci60_mod <- NCI60_4arrays
rownames(nci60_mod$hgu95) <- paste(rownames(nci60_mod$hgu95), "s1", sep="_")
mcoin_mod <- mcia(nci60_mod)
# without specifing
selectVar(mcoin_mod, a1.lim=c(2, Inf), a2.lim=c(-Inf, Inf))
# specifing the sepID.data and sepID.sep
selectVar(mcoin_mod, a1.lim=c(2, Inf), a2.lim=c(-Inf, Inf), sepID.data=4, sepID.sep="_")


}