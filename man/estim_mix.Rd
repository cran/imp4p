\name{estim.mix}
\alias{estim.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Estimation of a mixture model of MCAR and MNAR values in each column of a data matrix.
}
\description{
This function allows estimating a mixture model of MCAR and MNAR values in each column of data sets similar to the ones which can be studied in MS-based quantitative proteomics. Such data matrices contain intensity values of identified peptides.
}
\usage{
estim.mix(tab, tab.imp, conditions, x.min=20, x.max=30, x.step.mod=300, 
x.step.pi=300, nb.rei=100, method=2, gridsize=300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{tab.imp}{A matrix where the missing values of \code{tab} have been imputed under the assumption that they are all MCAR. For instance, such a matrix can be obtained by using the function \code{\link{impute.slsa}} of this package.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{x.min}{The lower bound of the interval used for estimating the cumulative distribution functions of the mixing model in each column.}
  \item{x.max}{The upper bound of the interval used for estimating the cumulative distribution functions of the mixing model in each column.}
  \item{x.step.mod}{The number of points in the intervals used for estimating the cumulative distribution functions of the mixing model in each column.}
  \item{x.step.pi}{The number of points in the intervals used for estimating the proportion of MCAR values in each column.}
  \item{nb.rei}{The number of initializations of the minimization algorithm used to estimate the proportion of MCAR values (see Details).}
  \item{method}{A numeric value indicating the method to use for estimating the proportion of MCAR values (see Details).}
  \item{gridsize}{A numeric value indicating the number of possible choices used for estimating the proportion of MCAR values with the method of Patra and Sen (2015) (see Details).}
}
\details{
This function aims to estimate the following mixture model in each column:

\code{F_tot(x)=pi_na*F_na(x)+(1-pi_na)*F_obs(x)}

\code{F_na(x)=pi_mcar*F_tot(x)+(1-pi_mcar)*F_mnar(x)}

where \code{pi_na} is the proportion of missing values, \code{pi_mcar} is the proportion of MCAR values, \code{F_tot} is the cumulative distribution function (cdf) of the complete values, \code{F_na} is the cdf of the missing values, \code{F_obs} is the cdf of the observed values, and \code{F_mnar} is the cdf of the MNAR values. 

To estimate this model, a first step consists to compute a rough estimate of \code{F_na(x)}  by assuming that all missing values are MCAR (thanks to the argument \code{tab.imp}). This rough estimate is noted \code{hat{F}_na(x)}. 

In a second step, the proportion of MCAR values is estimated. To do so, the ratio

\code{hat{pi}(x)=(1-hat{F}_na(x))/(1-hat{F}_tot(x))} 

is computed for different \code{x}, where 

\code{hat{F}_tot(x)=pi_na*hat{F}_na(x)+(1-pi_na)*hat{F}_obs(x)}

with \code{hat{F}_obs(x)} the empirical cdf of the observed values.

Next, the following minimization is performed:

\code{min_{1>k>0,a>0,d>0}} 

\code{sum_x (1/(s(x))^2)*[hat{pi}(x)-k-(1-k)*exp(-a*[(x-lower_x)]^d)/(1-hat{F}_tot(x))]^2}

where \code{(s(x))^2} is an estimate of the asymptotic variance of \code{hat{pi}(x)}, \code{lower_x} is an estimate of the minimum of the complete values and \code{upper_x} is the minimum between an estimate of the maximum of the missing values and the maximum of the observed values. To perform this minimization, the function \code{optim} with method "L-BFGS-B" is used. Because it is function of its initialization, it is possible to reinitialize a number of times the minimisation algorithm with the argument \code{nb.rei}: the parameters leading to the lowest minimum are next kept. 

Once \code{k}, \code{a} and \code{d} are estimated, one can use several methods to estimate \code{pi_mcar}: it is estimated 

by \code{k} if \code{method=1}, 

by 

\code{k+(1-k)*exp(-a*[(max(x)-lower_x)]^d)/(1-hat{F}_tot(max(x)))} 

if \code{method=2}, 

by estimating a decreasing trend with the PAVA algorithm on 

\code{k+(1-k)*exp(-a*[(x-lower_x)]^d)/(1-hat{F}_tot(x))} 

and keeping the righmost value of this trend if \code{method=3}, 

by estimating a decreasing trend with the PAV algorithm on \code{hat{pi}(x)} and keeping the righmost value of this trend if \code{method=4}, 

by using the histogram of \code{hat{pi}(x)} if \code{method=5}, 

by using the method of Patra and Sen (2015) adapted to our problematic if \code{method=6}. 
}
\value{
A list composed of:
  \item{abs.pi}{A numeric matrix containing the intervals used for estimating the ratio 
  
  \code{(1-F_na(x))/(1-F_tot(x))} 
  
  in each column.}
  \item{pi.init}{A numeric matrix containing the estimated ratios 
  
  \code{(1-F_na(x))/(1-F_tot(x))} 
  
  where \code{x} belongs to \code{abs.pi[,j]} for each sample \code{j}.}
  \item{var.pi.init}{A numeric matrix containing the estimated asymptotic variances of \code{pi.init}.}
  \item{abs.mod}{A numeric vector containing the interval used for estimating the mixture models in each column.}
  \item{pi.na}{A numeric vector containing the proportions of missing values in each column.}
  \item{F.na}{A numeric matrix containing the estimated cumulative distribution functions of missing values in each column on the interval \code{abs.mod}.}
  \item{F.tot}{A numeric matrix containing the estimated cumulative distribution functions of complete values in each column on the interval \code{abs.mod}.}
  \item{F.obs}{A numeric matrix containing the estimated cumulative distribution functions of observed values in each column on the interval \code{abs.mod}.}
  \item{F.mnar}{A numeric matrix containing the estimated cumulative distribution functions of MNAR values in each column on the interval \code{abs.mod}.}
  \item{pi.mcar}{A numeric vector containing the estimations of the proportion of MCAR values in each column.}
}
\references{
Patra, R. K., & Sen, B. (2015). Estimation of a two component mixture model with applications to multiple testing. Journal of the Royal Statistical Society: Series B (Statistical Methodology).
}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{impute.slsa}}
}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600,pi.mcar=0.2,para=10,nb.cond=2,nb.repbio=3,
nb.sample=5,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2);

#Deleting rows without any observed value in a condition
result=delete.na.rows(tab=res.sim$dat.obs, tab.c=res.sim$dat.comp, conditions=res.sim$conditions,
list.MCAR=res.sim$list.MCAR);

#Imputation of missing values with the slsa algorithm
dat.slsa=impute.slsa(result$tab.mod,conditions=res.sim$conditions,repbio=res.sim$repbio,
reptech=as.factor(1:length(res.sim$conditions)));

#Estimation of the mixture model
res=estim.mix(tab=result$tab.mod, tab.imp=dat.slsa, conditions=res.sim$condition);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
