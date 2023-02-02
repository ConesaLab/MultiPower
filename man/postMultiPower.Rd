\name{postMultiPower}
\alias{postMultiPower}
\title{
Power study for different sample sizes.
}
\description{
When the optimal sample size estimated by MultiPower exceeds the available budget, an alterative solution is to decrease such sample size at the cost of removing features
  with low effect size (magnitude of change). This function answer the question of to which effect size users must restrict themselves for a given maximum sample size. 
    Pilot multi-omic data sets are needed to perform this analysis. 
}
\usage{
postMultiPower(optResults, max.size = 5, omicCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
  \item{optResults}{
R object containing the results of \code{MultiPower} function.
}
  \item{max.size}{
Maximum sample size allowed by the user. It will be used to determine the effect size that can be detected (by default, 5).
}
  \item{omicCol}{
The color that will be used to plot each omic. It must be a vector with length equal to the number of omics. If it is NULL (default), default colors are used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a list with the following elements:
\item{SampleSize }{Matrix containing the optimal sample size for each omic data type (in columns) and for different cutoffs of Cohen’s d (in rows).}
\item{Power }{Matrix containing the statistical power at the optimal sample size for each omic data type (in columns) and for different cutoffs of Cohen’s d (in rows).}
\item{NumFeat }{Matrix containing the number of remaining features for each omic data type (in columns) and for different cutoffs of Cohen’s d (in rows).}
\item{d }{Values of Cohen’s d used as cutoffs to remove low effect size features.}
In addition, it also returns two plots that summarize these results: first, the number of replicates (sample size) for each of the tested Cohen's d values; 
second, the statistical power for each omic at the previously obtained optimal sample size and at the current sample size set by the user.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
