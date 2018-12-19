\name{MultiPower}
\alias{MultiPower}
\title{
Optimal sample size estimation and power study.
}
\description{
MultiPower estimates the optimal sample size for a multi-omic experiment when pilot multi-omic data sets are available to estimate the parameters required to compute power. An optimization problem is solved in order to achieve the desired power while minimizing the cost of the experiment.
}
\usage{
MultiPower(data, groups, type, d0 = 0.8, p1 = 0.2, omicPower = 0.6, averagePower = 0.85, dispPerc = 75,
fdr = 0.05, alpha = 0.05, cost = 1, equalSize = TRUE, max.size = 200, omicCol = NULL, powerPlots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
List with as many elements as omic data types. The names of the omics should be the names of the list. Each element in this list must be a raw count data matrix, and in this case MultiPower will take into account the library sizes to estimate power; a normally distributed data matrix which must have been already pre-processed and normalized; or a binary data matrix (with 0/1 or TRUE/FALSE values). In any case, for each one of these matrices, rows must correspond to omic features (genes, methylation sites, ChIP-seq peaks, etc.) and columns to observations (biological samples, patients, etc.).
}
  \item{groups}{
List with as many elements as omic data types. The names of the omics should be the names of the list. Each element in this list must be a vector with length equal to the number of observations for that omic in data argument. Each element of this vector must indicate the experimental group where each observation belong. Only two experimental groups are allowed.
}
  \item{type}{
Vector with length equal to the number of omic data types. Each element of this vector must be a 1, 2 or 3 to indicate whether the omic data are count data (1), continuous data approximately following a normal distribution (2) or binary data (3).
}
  \item{d0}{
Initial Cohen's d values. It can be either a single numeric value if the same Cohen's d is required for all the omics, or a vector with length equal to the number of omics containing the targeted Cohen’s d value for each omic. The Cohen's d is defined as delta/sigma, where delta is the absolute difference of means for each experimental group and sigma is the pooled standard deviation of both groups. The Cohen's d is not dependent on the scale of the data as it happens with  value, so it is fair to set the same value for all the omics. MultiPower computes the Cohen's d value for all the omic features and uses it to estimate the set M1 of differentially expressed (DE) features (those with d>d0), which are called pseudo-DE features. By default, d0 = 0.8. Note that for binary data, d0 is the Cohen's h value.
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
  \item{dispPerc}{
Dispersion percentile to be used to estimate power. By default, 75 (75th percentile).
}
  \item{fdr}{
False Discovery Rate level to be used. It is the significance level after multiple testing correction. By default, fdr = 0.05. If no multiple testing correction is to be applied, this argument must be set to NULL and then alpha argument is required.
}
  \item{alpha}{
Significance level to be used only when multiple testing correction is not to be applied (fdr = NULL). By default, alpha = 0.05.
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
  \item{powerPlots}{
If TRUE (default), power plots will be generated.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
When applying \code{MultiPower}, the result is a list containing the following elements:
\item{parameters }{List with as many elements as omic data types. For each omic, each element of the list is another list containing the different parameters used to compute power, either estimated from the pilot data or provided by the user: type, allDispersions, dispersion, p1, d, delta, minFC, meanCounts, m, m1, and alld.}
\item{optimalSampleSize }{List containing the following elements: n0 (sample size to achieve the minimum omic power, omicPower, for each omic), n (optimal sample size), finalPower (power at the optimal sample size for each omic), fdr (see fdr argument), omicPower (see omicPower argument), averagePower (see averagePower argument), and cost (see cost argument).}
\item{summary }{Table summarizing MultiPower results. The columns are: the names of the omic data sets (omic), the omic data type (type), the number of omic features for each omic (numFeat), the expected proportion of DE features set by the user (DEperc), the initial Cohen’s d set by the user (CohenD), the estimated delta value (delta), which is the absolute difference of group means for normal data and the proportions per group in binary data, the estimated fold-change in count data (minFC), the estimated mean of counts for count data (meanCounts), the estimated dispersion (dispersion), the FDR value (FDR), the minimum power to be achieve for each omic (minPower), the average power to be achieved in the multi-omic experiment (averPower), the cost per omic (cost), the minimum sample size needed for each omic to achieve minPower (minSampleSize), the optimal sample size (optSampleSize), and the power at this optimal sample size (power).}
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
