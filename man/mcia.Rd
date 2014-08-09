\name{mcia}
\alias{mcia}
\alias{plot.mcia}
\title{multiple co-inertia analysis}

\description{
    The main function in omicade4.
    Performing multiple co-inertia analysis on a \code{list} of \code{data.frame}s 
    or \code{matrix}
}

\usage{
mcia(df.list, cia.nf = 2, cia.scan = FALSE, nsc = T, svd = TRUE)
\method{plot}{mcia}(x, axes = 1:2, 
  sample.lab = TRUE, sample.legend = FALSE, sample.color = 1, 
  phenovec = NULL, df.color = 1,
  df.pch = NA, gene.nlab = 0, ...)
}

\arguments{
  \item{df.list}{
    A list of \code{data.frame}s, \code{matrix} or \code{ExpressionSet} is going to be analyzed, 
    the column number must be the same and mapped across all \code{data.frame}/\code{matrix}
}
  \item{cia.nf}{
    An integer indicating the number of kept axes
}
  \item{cia.scan}{
    A logical indicating whether the co-inertia analysis
    eigenvalue (scree) plot should be shown so that the number of axes, 
    (\code{cia.nf}) can be selected interactively. Default value is FALSE.
}
  \item{nsc}{
    A logical indicating whether multiple co-inertia analysis should be
    performed using multiple non-symmetric correspondence analyses 
    \code{\link[ade4:dudi.nsc]{dudi.nsc}}. The default =TRUE is highly 
    recommended. If FALSE, COA \code{\link[ade4:dudi.coa]{dudi.coa}} 
    will be performed on the first \code{data.frame}, and row weighted 
    COA \code{\link[made4:dudi.rwcoa]{dudi.rwcoa}} will be performed on 
    the rest ones using the row weights from the first one.
}
  \item{svd}{
    A logical indicates which function should be used to perform singular value decomposition.
}
  \item{sample.lab}{
    A logical indicating if the samples should be labelled, the default is
    TRUE.
}
  \item{sample.color}{
    Defining colours of samples for plotting sample space, the length of this 
    argument should be either one (uniform color) or the same with the 
    column number of \code{data.frame} in \code{df.list}. When sample
    size is big, \code{phenovec} should be used to distinguish 
    subgroups (See examples).
}
  \item{sample.legend}{
    A logical indicating if the legend for sample space should be drawn
    or a charater string passed to draw self-defined legend (see example).
}
  \item{df.color}{
    Defining the colours for plotting variables (genes) from different \code{data.frame}.
    The length of this argument should be either one (all datasets use the same
    colour) or the same number of datasets (each dataset has a specified colour, the
    repetitive use of colour code is allowed.)
}
  \item{df.pch}{
    Defining the \code{pch} for plotting variable (gene) space. The default is NA, the function
    will distinguish datasets by default. Otherwise, the length of this argument 
    should be either one (all datasets use the same \code{pch}) or the same number 
    of datasets (each dataset has a specified \code{pch}).
}
  \item{phenovec}{
    A factor for plotting sample space, phenovec could be 
    used to distinguish individuals in the \code{data.frame}s (see examples).
}
  \item{x}{
    An object of class \code{mcia}
}
  \item{axes}{
    A vector of integer in length 2 to indicate the axes are going to be plotted.
    The default are first two axes.
}
  \item{gene.nlab}{
    An integer indicating how many top weighted genes on each axis should be labelled
}
  \item{\dots}{
    Other arguments
}
}
\details{
    The column number of \code{data.frame} in the \code{df.list} must be the same, 
    and the same column
    from different \code{data.frame} should be matchable. For example, Microarray
    profiling for the same set of cell lines, patients and etc. 

    \code{mcia} calls \code{\link[ade4:dudi.nsc]{dudi.nsc}}, 
    \code{\link[ade4:ktab]{ktab}} and \code{\link[ade4:mcoa]{mcoa}} in ADE4
    packages.

    \bold{Plotting and visualizing \code{mcia} results}

    Two functions could be used to visualize the result of \code{mcia}:
    The first is \code{plot.mcia},
    which results in four plots. Top left represents the sample space. Individuals
    from the same column of different \code{data.frame}s are linked by edges.
    Different platforms are distinguished by the shape of points. 
    Top right shows the variable space, datasets are marked by different colours.
    Bottom left represents the eigenvalue scree plot. 
    The pseudo-eigenvalue space of all \code{data.frame}s are visualized in the bottom right panel. 
    The second function is \code{plotVar.mcia}, which could be used to
    plot the variable space for different datasets as well as finding and visualizing the
    variables (genes) across datasets.

    \bold{Other methods}

    \code{selectVar.mcia}: selecting variables (genes) according to the their coordinates.    
}

\value{
  \item{call }{
    the function called
}
  \item{mcoa }{
    The results returned by \code{\link[ade4:mcoa]{mcoa}}
}
  \item{coa }{
    The results returned by separate analysis (applying \code{\link[ade4:dudi.nsc]{dudi.nsc}}
    or \code{\link[ade4:dudi.coa]{dudi.coa}} on each \code{data.frame} separately)
}
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
  See Also as \code{\link[ade4:mcoa]{mcoa}}, \code{plotVar}, \code{plotVar}
}
\examples{

data(NCI60_4arrays)
mcoin <- mcia(NCI60_4arrays)
plot(mcoin, sample.lab=FALSE, df.col=4:7)

colcode <- sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\\."), "[", 1)

# "sample.color" should be used when sample number is small, this is not proper use here 
plot(mcoin, sample.lab=FALSE, sample.color=as.factor(colcode), sample.legend=TRUE)
# phenovec should be used to use color distinguishing subtypes/groups (when sample size is large)
plot(mcoin, sample.lab=FALSE, phenovec=as.factor(colcode), sample.legend=TRUE)
# or use self defined legend
leg <- 'legend("topleft", legend=unique(colcode), pch=20, col=1:9, border=FALSE, bty="n")'
plot(mcoin, sample.lab=FALSE, sample.legend=leg, sample.color=as.factor(colcode))

}


\keyword{mcia}
