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
False Discovery Rate level, that is significance level after multiple testing correction, to be used. By default, fdr = 0.05.
}
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
myomics = c("RNA-seq", "miRNA-seq",  "RRBS-seq")
areCounts = c(TRUE, TRUE, FALSE)
sdValue = c(2.5,3.7,0.7)  # pooled standard deviation
p1omics = c(0.4, 0.2, 0.1) # expected proportion of DE features
d0 = 1
minFC0 = c(1.5, 1.5, NA) # fold-change to be detected in DE features
meanCounts0 = c(20, 15, NA) # mean counts of DE features
M = c(12000, 500, 600000) # number of omic features
myparam = vector("list", length = 3)  # define parameters object
names(myparam) = myomics
for (i in 1:length(myparam)) {
  if (areCounts[i]) {
    disp = max((sdValue[i]^2 - meanCounts0[i]) / meanCounts0[i]^2, 0.0001)
  } else { disp = sdValue[i] }
  myparam[[i]] = list("counts" = areCounts[i], "allDispersions" = NULL,
"dispersion" = disp,"p1" = p1omics[i], "d" = d0,
"delta" = d0*sdValue[i], "minFC" = minFC0[i],
"meanCounts" = meanCounts0[i], "m" = M[i],
"m1" = round(M[i]*p1omics[i], 0), "alld" = NULL)
}
optimalSS = optimalRep(parameters = myparam, omicPower = 0.6,
averagePower = 0.8, fdr = 0.05, cost = 1,
equalSize = TRUE, max.size = 30)
optimalSS$n  # optimal sample size
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
