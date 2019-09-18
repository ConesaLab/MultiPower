\name{MultiCompaPlot}
\alias{MultiCompaPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary plots for multiple comparisons.
}
\description{
Plotting sample size and statistical power per omic for multiple pairwise comparisons.
}
\usage{
MultiCompaPlot(multiOutput, omicCol = NULL, equalSize, legendLoc = "bottomright")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiOutput}{
R object containing the results of \code{MultiGroupPower}.
}
  \item{omicCol}{
The color that will be used to plot each omic. It must be a vector with length equal to the number of omics. If it is NULL (default), default colors are used.
}
  \item{equalSize}{
If TRUE, the same sample size was required for all omics when applying \code{MultiGroupPower} to obtain the results in multiOutput argument. If FALSE, a different sample size was allowed for each omic data type.
}
  \item{legendLoc}{
Where to locate the legend inside the plots. By default, "bottomright". The options are "bottomright", "bottom",
"bottomleft", "left", "topleft", "top", "topright", "right" and "center".
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%
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
