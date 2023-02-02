\name{optimalRep}
\alias{optimalRep}
\title{
Estimation of optimal sample size.
}
\description{
Estimation of the optimal sample size when pilot multi-omic data sets are not available.
}
\usage{
optimalRep(parameters, omicPower = 0.6, averagePower = 0.85, fdr = 0.05, cost = 1, equalSize = TRUE, max.size = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
List with as many elements as omic data types. For each omic, each element of this list is another list containing the different parameters needed to compute power which, in this case, must be set by the user. See \code{MultiPower} for more details.
}
  \item{omicPower}{
The minimum power that must be achieved for each omic. It must be a vector with length equal to the number of omics. If it is a single number, this same number will be used for all the omics. By default, omicPower = 0.6.
}
  \item{averagePower}{
The minimum average power that must be globally achieved. By default, averagePower = 0.85.
}
  \item{fdr}{
False Discovery Rate level to be used. It is the significance level after multiple testing correction. By default, fdr = 0.05.}
  \item{cost}{
The cost to generate a replicate (a sample) for each omic. It must be a vector with length equal to the number of omics. If it is a single number, this same number will be used for all the omics. This argument will only be used when a different sample size per omic is allowed. By default, cost = 1 (which means that all the omics will be assumed to have the same cost).
}
  \item{equalSize}{
 If TRUE (default), the same optimal sample size will be estimated for all the omics. If FALSE, omics are allowed to have different sample sizes.
}
  \item{max.size}{
Maximum allowed sample size. By default, max.size = 30.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
Sonia Tarazona; David Gómez-Cabrero
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MultiPower}}
}
\examples{
optimalSS = optimalRep(parameters = myparam, omicPower = 0.6,
averagePower = 0.8, fdr = 0.05, 
cost = 1, equalSize = TRUE, max.size = 30)
optimalSS$n  # optimal sample size
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
