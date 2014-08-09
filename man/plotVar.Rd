\name{plotVar}
\alias{plotVar}
\title{
  Plot variable (gene) spaces of result from MCIA or CIA
}
\description{
  The user level function for plotting variable space of \code{\link{mcia}} 
  or \code{\link[made4:cia]{cia}}, which could be used to visualize selected 
  variables (genes) across datasets. It calls \code{plotVar.cia} or \code{plotVar.mcia}.
}
\usage{
  plotVar(x, var = NA, axes = 1:2, 
          var.col = "red", var.lab = FALSE, bg.var.col = "gray", 
          nlab = 0, sepID.data=NULL, sepID.sep="_", ...)
}
\arguments{
  \item{x}{
    An object of class \code{cia} or \code{mcia}
}
  \item{var}{
    A character vector defining the variables (genes) are going to be labelled
    and coloured. The default NA means no variables (genes) selected.
}
  \item{axes}{
    An integer vector in length 2 indicating which axes are going to be plotted. 
    Default are the first two axes.
}
  \item{var.col}{
    The colour of selected variables (genes), the length of this argument should be 
    either 1 (uniform colour) or the length of \code{var} (each \code{var}
    has a specified colour).
}
  \item{var.lab}{
    A logical indicating if the variables (genes) selected should be labelled, the 
    default is FALSE
}
  \item{bg.var.col}{
    Colour code for unselected variables (genes) in all datasets.
}
  \item{nlab}{
    An integer indicating how many top weighted genes on each axis should be labelled.
}
  \item{sepID.data}{
    This argument enables a more generalized mapping of identifiers in different datasets.
    For example, if there is a PTM (post-transcriptional modification) dataset in one of 
    the \code{data.frame}s, the corresponding protein could 
    be detected with setting this argument. For more details, see "details" section.
}
  \item{sepID.sep}{
    Used to help determine the separator of variables (genes) in the sepID.data. 
    For more details, see "details" section.
}
  \item{\dots}{
    Other arguments
}
}
\details{
  For the sepID.data, a typical example is the post-transcriptional modification (PTM) data.
  The name of variables (genes) have a general form like 
  "proteinName_modificationSite". The \code{sepID.data} specifies the IDs from dataset
  that should be separated, \code{sepID.sep} specifies the separator of protein name
  and modification site. This is used to determine the same proteins/genes
  across different datasets.
}
\value{
  If \code{var} is not NA, a data frame is returned, with rows for variables (genes) of
  interest and columns of logical values indicating which dataset contains which
  variables (genes). 
}
%\references{
%  Meng C, Kuster B, Culhane AC and Gholami AM. A multivariate 
%  approach to the integration of multi-omics datasets. 
%  BMC bioinformatics (under review)
%}
\author{
  Chen Meng
}


\seealso{
  See Also as \code{\link{plotVar.cia}}, \code{\link{plotVar.mcia}}
}

\examples{

data(NCI60_4arrays)
mcoin <- mcia(NCI60_4arrays)
plotVar(mcoin, var=c("S100B", "S100A1"), var.lab=TRUE)

# an example for the usage of sepID.data and sepID.sep
nci60_mod <- NCI60_4arrays
rownames(nci60_mod$hgu95) <- paste(rownames(nci60_mod$hgu95), "s1", sep="_")
mcoin_mod <- mcia(nci60_mod)
# without specifying
plotVar(mcoin_mod, var=c("S100B", "S100A1"), var.lab=TRUE)
# specifying the sepID.data and sepID.sep
plotVar(mcoin_mod, var=c("S100B", "S100A1"), var.lab=TRUE, sepID.data=4, sepID.sep="_")

}
