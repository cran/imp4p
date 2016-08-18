\name{impute.mix}
\alias{impute.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation under an assumption of a mixture of MCAR and MNAR values. 
}
\description{
This function allows imputing data sets with a MCAR-devoted algorithm and a MNAR-devoted algorithm using the probabilities that missing values are MCAR.
}
\usage{
impute.mix(tab, prob.MCAR, conditions, repbio, reptech, method="slsa", nknn=15,
weight=1, selec="all", ind.comp=1, progress.bar=TRUE, q=0.95)
 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
\item{prob.MCAR}{A matrix of probabilities that each missing value is MCAR. For instance such a matrix can be obtained from the function \code{\link{prob.mcar.tab}} of this package.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{repbio}{A vector of factors indicating the biological replicate to which each column belongs.}
  \item{reptech}{A vector of factors indicating the technical replicate to which each column belongs.}
  \item{method}{The method used for imputing MCAR data. If \code{methodi="slsa"} (default), then the SLSA algorithm is used, else the MLE algorithm is used.}
  \item{nknn}{The number of nearest neighbours used in the SLSA algorithm (see \code{\link{impute.slsa}}).}
  \item{weight}{The way of weighting in the algorithm (see \code{\link{impute.slsa}}).}
  \item{selec}{A parameter to select a part of the dataset to find nearest neighbours between rows. This can be useful for big data sets (see \code{\link{impute.slsa}}).}
  \item{ind.comp}{If \code{ind.comp=1}, only nearest neighbours without missing values are selected to fit linear models (see \code{\link{impute.slsa}}). Else, they can contain missing values.}
  \item{progress.bar}{If \code{TRUE}, a progress bar is displayed.}
    \item{q}{A quantile value (see \code{\link{impute.igcda}}).}
}
\details{
The missing values for which \code{prob.MCAR} is superior to 0.5 are imputed with either the function \code{\link{impute.slsa}} or the MLE algorithm (function \code{impute.wrapper.MLE} of the R package imputeLCMD). The other missing values are considered MNAR and imputed with \code{\link{impute.igcda}}. 
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

#Deleting rows without any observed value in a condition
result=delete.na.rows(tab=res.sim$dat.obs, tab.c=res.sim$dat.comp, conditions=res.sim$conditions,
list.MCAR=res.sim$list.MCAR);

#Imputation of missing values with the slsa algorithm
dat.slsa=impute.slsa(result$tab.mod,conditions=res.sim$condition,repbio=res.sim$repbio,
reptech=as.factor(1:length(res.sim$condition)));

#Estimation of the mixture model
res=estim.mix(tab=result$tab.mod, tab.imp=dat.slsa, conditions=res.sim$condition);

#Computing probabilities to be MCAR
born=estim.bound(tab=result$tab.mod,conditions=res.sim$condition);
proba=prob.mcar.tab(born$tab.lower,born$tab.upper,res);

#Imputation under the assumption of MCAR and MNAR values
tabi=impute.mix(tab=result$tab.mod, prob.MCAR=proba, conditions=res.sim$conditions,
repbio=res.sim$repbio, reptech=as.factor(1:length(res.sim$condition)), method="slsa",
nknn=15, weight=1, selec="all", ind.comp=1, progress.bar=TRUE);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line