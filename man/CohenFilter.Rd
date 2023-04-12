\name{CohenFilter}
\alias{CohenFilter}
\title{
Removal of low effect features.
}
\description{
Removal of omic features with Cohen’s d (or h) below the threshold set by the user.
}
\usage{
CohenFilter(data, d, parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
Original data given to MultiPower function.
}
  \item{d}{
Cutoff value for Cohen’s d (or h). A scalar (and then the same cutoff is applied to all the omics) or a vector with as many values as omic data types in the data object.
}
  \item{parameters}{
List with as many elements as omic data types estimated by MultiPower function. They can be retrieved from MultiPower output as output@parameters.
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
%%
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
