\name{MultiPower}
\alias{MultiPower}
\title{
Optimal sample size estimation and power study.
}
\description{
MultiPower estimates the optimal sample size for a multi-omic experiment when pilot multi-omic data sets are available to estimate the parameters required to compute power. An optimization problem is solved in order to achieve the desired power while minimizing the cost of the experiment.
}
\usage{
MultiPower(data, groups, counts, d0 = 0.8, p1 = 0.2, omicPower = 0.6, averagePower = 0.85,
fdr = 0.05, cost = 1, equalSize = TRUE, max.size = 200, omicCol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
List with as many elements as omic data types. The names of the omics should be the names of the list. Each element in this list must be the already pre-processed and normalized data matrix of an omic data type, where rows correspond to omic features (genes, methylation sites, ChIP-seq peaks, etc.) and columns to observations (biological samples, patients, etc.).
}
  \item{groups}{
List with as many elements as omic data types. The names of the omics should be the names of the list. Each element in this list must be a vector with length equal to the number of observations for that omic in data argument. Each element of this vector must indicate the experimental group where each observation belong. Only two experimental groups are allowed.
}
  \item{counts}{
Vector with length equal to the number of omic data types. Each element of this vector must indicate whether the omic data are count data (TRUE) or continuous data approximately following a normal distribution (FALSE).
}
  \item{d0}{
Initial Cohen’s d value. The Cohen’s d is defined as delta/sigma, where delta is the absolute difference of means for each experimental group and sigma is the pooled standard deviation of both groups. The Cohen’s d is not dependent on the scale of the data as it happens with delta value, so it is fair to set the same value for all the omics. By default, d0 = 0.8.
}
  \item{p1}{
The expected proportion of DE features. It must be a vector with length equal to the number of omics. If it is a single number, this same number will be used for all the omics. The values for p1 must be set by the user, who can decide to estimate them from pilot data by performing a differential expression analysis (this option is not provided by MultiPower). By default, p1 = 0.2.
}
  \item{omicPower}{
The minimum power that must be achieved for each omic. It must be a vector with length equal to the number of omics. If it is a single number, this same number will be used for all the omics. By default, omicPower = 0.6.
}
  \item{averagePower}{
The minimum average power that must be globally achieved. By default, averagePower = 0.85.
}
  \item{fdr}{
False Discovery Rate level, that is significance level after multiple testing correction, to be used. By default, fdr = 0.05.
}
  \item{cost}{
The cost to generate a replicate (a sample) for each omic. It must be a vector with length equal to the number of omics. If it is a single number, this same number will be used for all the omics. This argument will only be used when a different sample size per omic is allowed. By default, cost = 1 (which means that all the omics will be assumed to have the same cost).
}
  \item{equalSize}{
If TRUE (default), the same optimal sample size will be estimated for all the omics. If FALSE, omics are allowed to have different sample sizes.
}
  \item{max.size}{
Maximum allowed sample size. By default, max.size = 200.
}
  \item{omicCol}{
The color that will be used to plot each omic. It must be a vector with length equal to the number of omics. If it is NULL (default), default colors are used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
When applying \code{MultiPower}, the result is a list containing the following elements:
\item{parameters }{List with as many elements as omic data types. For each omic, each element of the list is another list containing the different parameters used to compute power, either estimated from the pilot data or provided by the user: counts, allDispersions, dispersion, p1, d, delta, minFC, meanCounts, m, m1, and alld.}
\item{optimalSampleSize }{List containing the following elements: n0 (sample size to achieve the minimum omic power, omicPower, for each omic), n (optimal sample size), finalPower (power at the optimal sample size for each omic), fdr (see fdr argument), omicPower (see omicPower argument), averagePower (see averagePower argument), and cost (see cost argument).}
\item{summary }{Table summarizing MultiPower results. The columns are: the omic data types (omic), whether the omic data are count data or not (counts), the number of omic features for each omic (numFeat), the expected proportion of DE features set by the user (DEperc), the initial Cohen’s d set by the user (CohenD), the estimated delta value, i.e. the absolute difference of group means (delta), the estimated fold-change in count data (minFC), the estimated mean of counts for count data (meanCounts), the estimated dispersion (dispersion), the FDR value (FDR), the minimum power to be achieve for each omic (minPower), the average power to be achieved in the multi-omic experiment (averPower), the cost per omic (cost), the minimum sample size needed for each omic to achieve minPower (minSampleSize), the optimal sample size (optSampleSize), and the power at this optimal sample size (power).}
\item{data2plot}{Data generated to create the power plots that are also returned by the function.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sonia Tarazona; David Gómez-Cabrero
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
