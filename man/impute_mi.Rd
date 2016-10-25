\name{impute.mi}
\alias{impute.mi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of data sets containing peptide intensities with a multiple imputation strategy.
}
\description{
This function allows imputing data sets containing peptide intensities with a multiple imputation strategy.
}
\usage{
impute.mi(tab, conditions, repbio=NULL, reptech=NULL, nb.iter=3, nknn=15, selec=600, 
siz=500, weight=1, ind.comp=1, progress.bar=TRUE, x.min=20, x.max=30, x.step.mod=300, 
x.step.pi=300, nb.rei=100, method=4, gridsize=300, q=0.95, q.min=0, q.norm=3, 
eps=2, methodi="slsa");
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing only numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{repbio}{A vector of factors indicating the biological replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{reptech}{A vector of factors indicating the technical replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{nb.iter}{The number of iterations used for the multiple imputation method (see \code{\link{mi.mix}}).}
  \item{methodi}{The method used for imputing data. If \code{methodi="mle"}, then the MLE algorithm is used (function \code{impute.wrapper.MLE} of the R package imputeLCMD), else the SLSA algorithm is used (default). (see \code{\link{mi.mix}})}
  \item{nknn}{The number of nearest neighbours used in the SLSA algorithm (see \code{\link{impute.slsa}}).}
  \item{selec}{A parameter to select a part of the dataset to find nearest neighbours between rows. This can be useful for big data sets (see \code{\link{impute.slsa}}).}
    \item{siz}{A parameter to select a part of the dataset to perform imputations with the SLSA algorithm or the MLE algorithm. This can be useful for big data sets (see \code{\link{mi.mix}}).}
  \item{weight}{The way of weighting in the algorithm (see \code{\link{impute.slsa}}).}
  \item{ind.comp}{If \code{ind.comp=1}, only nearest neighbours without missing values are selected to fit linear models (see \code{\link{impute.slsa}}). Else, they can contain missing values.}
  \item{progress.bar}{If \code{TRUE}, a progress bar is displayed.}
  \item{x.min}{The lower bound of the interval used for estimating the cumulative distribution functions of the mixing model in each column (see \code{\link{estim.mix}}).}
  \item{x.max}{The upper bound of the interval used for estimating the cumulative distribution functions of the mixing model in each column (see \code{\link{estim.mix}}).}
  \item{x.step.mod}{The number of points in the intervals used for estimating the cumulative distribution functions of the mixing model in each column (see \code{\link{estim.mix}}).}
  \item{x.step.pi}{The number of points in the intervals used for estimating the proportion of MCAR values in each column (see \code{\link{estim.mix}}).}
  \item{nb.rei}{The number of initializations of the minimization algorithm used to estimate the proportion of MCAR values (see Details) (see \code{\link{estim.mix}}).}
  \item{method}{A numeric value indicating the method to use for estimating the proportion of MCAR values (see \code{\link{estim.mix}}).}
  \item{gridsize}{A numeric value indicating the number of possible choices used for estimating the proportion of MCAR values with the method of Patra and Sen (2015) (see \code{\link{estim.mix}}).}
  \item{q}{A quantile value (see \code{\link{impute.igcda}}).}
  \item{q.min}{A quantile value of the observed values allowing defining the maximal value which can be generated. Default is 0 (the maximal value is the minimum of observed values minus \code{eps}) (see \code{\link{impute.pa}}).}
  \item{q.norm}{A quantile value of a normal distribution allowing defining the minimal value which can be generated. Default is 3 (the minimal value is the maximal value minus qn*median(sd(observed values)) where sd is the standard deviation of a row in a condition) (see \code{\link{impute.pa}}).}
  \item{eps}{A value allowing defining the maximal value which can be generated. Default is 2 (see \code{\link{impute.pa}}).}
}
\details{
First, a mixture model of MCAR and MNAR values is estimated in each column of \code{tab}. This model is used to estimate probabilities that each missing value is MCAR. Then, these probabilities are used to perform a multiple imputation strategy (see \code{\link{mi.mix}}). Rows with no value in a condition are imputed using the \code{\link{impute.pa}} function.
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

#Imputation of the dataset noting the conditions to which the samples belong.
result=impute.mi(tab=res.sim$dat.obs, conditions=res.sim$conditions);

#Imputation of the dataset noting the conditions to which the samples belong 
#and also their biological replicates.
result=impute.mi(tab=res.sim$dat.obs, conditions=res.sim$conditions, repbio=res.sim$repbio);

#For large data sets, the imputation can be accelerated thanks to the selec parameter 
#and the siz parameter (see impute.slsa and mi.mix)
#but it may result in a less accurate data imputation. Note that selec has to be greater than siz.
#
#Here, nb.iter is fixed to 3
result1=impute.mi(tab=res.sim$dat.obs, conditions=res.sim$conditions, progress.bar=TRUE, 
selec=400, siz=300, nb.iter=3);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
