\name{fuzzyROC}
\alias{fuzzyROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to produce fuzzyROC object for plotting
}
\description{
Function to produce fuzzyROC object for plotting from raw labeled biomarker data
}
\usage{
fuzzyROC(XX,YY,grayMax = .2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XX}{Numerical vector of biomarker values}
  \item{YY}{Binary vector of class labels}
  \item{grayMax}{Maximum size of gray area (aka indeterminacy zone), as a proportion of total the length of XX. This quantity is referred to as \eqn{ \gamma } in the companion article.}
}
\details{
The function assumes low levels of the marker predict class 0, i.e. the disease class.
Details and notation in the companion article.
}
\value{
\item{XX}{Returns Input XX}
\item{YY}{Returns Input YY}
\item{grayMax}{Returns Input grayMax}
\item{cMid}{Vector of midpoints \eqn{c_j} of candidate gray areas}
\item{cLstar}{Vector of lower limits \eqn{c_L^*}}
\item{cHstar}{Vector of higher limits \eqn{c_H^*}}
\item{grayStar}{Vector of proportions of points in the optimal gray areas.}
\item{F0star}{Vector of \eqn{F_0^*} values}
\item{F1star}{Vector of \eqn{F_1^*} values}
\item{F0dagger}{Vector of \eqn{F_0^{\dagger}} values}
\item{F1dagger}{Vector of \eqn{F_1^{\dagger}} values}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
