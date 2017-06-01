\name{fast_sim}
\alias{fast_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute similarity measures between a vector and each row of a matrix.
}
\description{
This function allows computing either the pairwise correlation between a vector and each row of a matrix if at least three side-by-side observed values, or the euclidean distance between side-by-side observed values if the length of the vector that is compared to each row of the matrix is less than 3. It is implemented thanks to the RCpp package.
}
\usage{
fast_sim(prot, mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prot}{
A numeric vector without missing values.
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
## The function is currently defined as
##function (prot, mat)
##{
##    .Call("imp4p_fast_sim", PACKAGE = "imp4p", prot, mat)
## }
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cpp wrapper function}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
