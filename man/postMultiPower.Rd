\name{postMultiPower}
\alias{postMultiPower}
\title{
Power study for different sample sizes.
}
\description{
When the optimal sample size estimated by MultiPower exceeds the available budget, an alterative solution is to decrease such sample size at the cost of detecting a larger effect size (magnitude of change). This function answer the question of to which effect size users must restrict themselves for a given maximum sample size. Pilot multi-omic data sets are needed to perform this analysis. Binary data is not accepted and Cohen's d must have been set to the same value for all the omics.
}
\usage{
postMultiPower(data, groups, optResults, max.size = 5, omicCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data provided to \code{MultiPower} function, that is, a List with as many elements as omic data types. The names of the omics should be the names of the list. Each element in this list must be a raw count data matrix, and in this case MultiPower will take into account the library sizes to estimate power; a normally distributed data matrix which must have been already pre-processed and normalized; or a binary data matrix (with 0/1 or TRUE/FALSE values). In any case, for each one of these matrices, rows must correspond to omic features (genes, methylation sites, ChIP-seq peaks, etc.) and columns to observations (biological samples, patients, etc.).}
  \item{groups}{
List with as many elements as omic data types. The names of the omics should be the names of the list. Each element in this list must be a vector with length equal to the number of observations for that omic in data argument. Each element of this vector must indicate the experimental group where each observation belong. Only two experimental groups are allowed.
}
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
\item{SampleSize }{Matrix containing the optimal sample size for each omic data type (in rows) and for different values of Cohen’s d (in columns).}
\item{Power }{Matrix containing the statistical power at the optimal sample size for each omic data type (in rows) and for different values of Cohen’s d (in columns).}
\item{d }{Values of Cohen’s d for which the optimal sample size has been estimated. These are the values represented in the columns of the above matrices.}
In addition, it also returns two plots that summarize these results: first, the number of replicates (sample size) for each of the tested Cohen's d values; second, the statistical power for each omic at the previously obtained optimal sample size and at the current sample size set by the user.
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
