\name{powerSummary}
\alias{powerSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary table with sample size optimization results.
}
\description{
This function returns a table that summarizes input parameters to estimate sample size and the output of the estimation including the optimal sample size for each omic. See Details for more information. It is used after applying \code{optimalRep}, when pilot multi-omic data sets are not available.
}
\usage{
powerSummary(parameters, optimalSampleSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
List with as many elements as omic data types. For each omic, each element of this list is another list containing the different parameters needed to compute power which, in this case, must be set by the user. See MultiPower results for more details.}
  \item{optimalSampleSize}{
R object containing the results of \code{optimalRep}.
}
}
\details{
The columns of the summary table are: the omic data names (omic), the omic data types (type), the number of omic features for each omic (numFeat), the expected proportion of DE features set by the user (DEperc), the initial Cohen’s d set by the user (CohenD), the estimated delta value (delta), i.e. the absolute difference of group means for normal data and the proportions per group in binary data, the estimated fold-change in count data (minFC), the estimated mean of counts for count data (meanCounts), the estimated dispersion (dispersion), the FDR value (FDR), the minimum power to be achieve for each omic (minPower), the average power to be achieved in the multi-omic experiment (averPower), the cost per omic (cost), the minimum sample size needed for each omic to achieve minPower (minSampleSize), the optimal sample size (optSampleSize), and the power at this optimal sample size (power).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{MultiPower}}
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
