\name{impute.pa}
\alias{impute.pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of peptides having no value in a biological condition (present in a condition / absent in another).
}
\description{
This function imputes missing values by small values. 
}
\usage{
impute.pa(tab, conditions, q.min=0, q.norm=3, eps=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{q.min}{A quantile value of the observed values allowing defining the maximal value which can be generated. This maximal value is defined by the quantile \code{q.min} of the observed values distribution minus \code{eps}. Default is 0 (the maximal value is the minimum of observed values minus \code{eps}).}
  \item{q.norm}{A quantile value of a normal distribution allowing defining the minimal value which can be generated. Default is 3 (the minimal value is the maximal value minus qn*median(sd(observed values)) where sd is the standard deviation of a row in a condition).}
  \item{eps}{A value allowing defining the maximal value which can be generated. This maximal value is defined by the quantile \code{q.min} of the observed values distribution minus \code{eps}. Default is 2.}
}
\details{
This function replaces the missing values of the rows by random draws from an uniform distribution between the defined minimal value and maximal value. The value of \code{eps} can be interpreted as a minimal fold-change value above which the present/absent peptides appear.
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
data.small.val=impute.pa(res.sim$dat.obs,res.sim$conditions);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
