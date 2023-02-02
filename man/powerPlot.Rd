\name{powerPlot}
\alias{powerPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power plots.
}
\description{
Plotting statistical power versus sample size per omic from optimal sample size results.
and statistical power versus Cohen's d cutoffs (in percentiles).
}
\usage{
powerPlot(parameters, optimalSampleSize, omicCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
List with as many elements as omic data types. For each omic, each element of this list is another list containing the different parameters needed to compute power. 
See \code{MultiPower} results for more details.
}
  \item{optimalSampleSize}{
R object containing the results of \code{optimalRep}.
}
  \item{omicCol}{
The color that will be used to plot each omic. It must be a vector with length equal to the number of omics. If it is NULL (default), default colors are used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List containing the values used for the plot or plots.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sonia Tarazona
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ plot }% use one of  RShowDoc("KEYWORDS")
%\keyword{ power }% __ONLY ONE__ keyword per line
}
