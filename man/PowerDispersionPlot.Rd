\name{PowerDispersionPlot}
\alias{PowerDispersionPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power vs Dispersion plot for a given sample size.
}
\description{
Plotting statistical power versus dispersion percentiles for given sample sizes. Pilot data is required for this plot.
}
\usage{
PowerDispersionPlot(n = 5, parameters, fdr= 0.05, alpha = 0.05, omicCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size available. It can be an integer if the same number of replicates is available for all the omic data types or a vector if different. By default, n = 5.
}
  \item{parameters}{
List with as many elements as omic data types. For each omic, each element of this list is another list containing the different parameters needed to compute power which, in this case, must be set by the user. See \code{MultiPower} results for more details.
}
  \item{fdr}{
False Discovery Rate level to be used. It is the significance level after multiple testing correction. By default, fdr = 0.05. If no multiple testing correction is to be applied, this argument must be set to NULL and then alpha argument is required.
}
  \item{alpha}{
Significance level to be used only when multiple testing correction is not to be applied (fdr = NULL). By default, alpha = 0.05.
}
  \item{omicCol}{
The color that will be used to plot each omic. It must be a vector with length equal to the number of omics. If it is NULL (default), default colors are used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix containing the values used for this plot.
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
