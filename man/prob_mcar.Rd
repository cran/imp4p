\name{prob.mcar}
\alias{prob.mcar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a vector of probabilities that missing values are MCAR.
}
\description{
This function returns a vector of probabilities that each missing value is MCAR from specified confidence intervals.
}
\usage{
prob.mcar(b.l,b.u,absc,pi.mcar,F.tot,F.na)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b.l}{A numeric vector of lower bounds for missing values.}
  \item{b.u}{A numeric vector of upper bounds for missing values.}
  \item{absc}{The interval on which is estimated the MCAR data mechanism.}
  \item{pi.mcar}{The estimated proportion of MCAR values.}
  \item{F.tot}{An estimation of the cumulative distribution function of the complete values on the interval \code{absc}.}
  \item{F.na}{An estimation of the cumulative distribution function of the missing values on the interval \code{absc}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector of estimated probabilities to be MCAR for missing values in the confidence intervals defined by \code{b.l} and \code{b.u}. The input arguments \code{absc}, \code{pi.mcar}, \code{F.tot} and \code{F.na} can be estimated thanks to the function \code{\link{estim.mix}}.
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

\seealso{
\code{\link{estim.mix}}
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

born=estim.bound(tab=result$tab.mod,conditions=res.sim$condition);

#Computing probabilities to be MCAR in the first column of result$tab.mod
proba=prob.mcar(b.l=born$tab.lower[,1],b.u=born$tab.upper[,1],absc=res$abs.mod,
pi.mcar=res$pi.mcar[1], F.tot=res$F.tot[,1], F.na=res$F.na[,1]);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
