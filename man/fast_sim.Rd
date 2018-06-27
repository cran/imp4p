\name{fast_sim}
\alias{fast_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute similarity measures between a vector and each row of a matrix.
}
\description{
This function allows computing either the pairwise correlation between a vector and each row of a matrix if each has at least 3 side-by-side observed values, or 1/(1+d) where d is the euclidean distance between side-by-side observed values if the vector contains less than 3 observed values. If the vector contains more than 3 observed values and a row of the matrix has less than 3 observed values, then it gives NA. It is implemented thanks to the RCpp package.
}
\usage{
fast_sim(prot, mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prot}{
A numeric vector containing numeric and missing values.
}
  \item{mat}{
A data matrix containing numeric and missing values.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector containing the values of the similarity measures between the \code{prot} vector and each row of the \code{mat} matrix.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{impute.slsa}}
%}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=20000,nb.miss=1000);

#Fast computation of similarities
fast_sim(res.sim$dat.obs[1,],res.sim$dat.obs);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cpp wrapper function}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
