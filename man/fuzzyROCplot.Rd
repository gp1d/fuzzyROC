\name{fuzzyROCplot}
\alias{fuzzyROCplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fuzzy ROC plotting
}
\description{
Fuzzy ROC plotting from the object created by the fuzzyROC function
}
\usage{fuzzyROCplot(fuzzyROCplot = function(fR,mainTitle="",set=TRUE,segm=FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fR}{Output of fuzzyROC function}
  \item{mainTitle}{Main title of plot}
  \item{set}{Logical. If TRUE the area between the bounds is colored gray.}
  \item{segm}{Logical. If TRUE segments between pairs of points corresponding to the same gray zone are connected by segments.}
}
\details{
}
\value{
A single plot.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Giovanni Parmigiani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## based on the example in the companion paper
## library(curatedOvarianData)
## data(GSE32063_eset)
## GeneName = "ZNF544"
## XX = exp ( exprs(GSE32063_eset)[GeneName,] )
## YY = 1 * ( pData(GSE32063_eset)[,"debulking"] == "optimal" )
## fROC = fuzzyROC(XX,YY,.1)
## fuzzyROCplot(fROC,mainTitle="Gamma = .10",set=FALSE,segm=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
