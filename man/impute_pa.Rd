\name{impute.pa}
\alias{impute.pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of peptides having no value in a biological condition.
}
\description{
This function imputes missing values by small values. 
}
\usage{
impute.pa(tab, q=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{q}{A quantile value.}
}
\details{
This function replaces the missing values of the rows by random draws from an uniform distribution between the minimum of observed values in the column and the \code{q}\% quantile of these observed values. 
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600,pi.mcar=0.2,para=10,nb.cond=2,nb.repbio=3,
nb.sample=5,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2);

#Imputation of the simulated data set with small values
data.small.val=impute.pa(res.sim$dat.obs);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
