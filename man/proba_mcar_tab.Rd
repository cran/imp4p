\name{prob.mcar.tab}
\alias{prob.mcar.tab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a matrix of probabilities that missing values are MCAR.
}
\description{
This function returns a matrix of probabilities that each missing value is MCAR from specified confidence intervals.
}
\usage{
prob.mcar.tab(tab.l,tab.u,res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab.l}{A numeric matrix of lower bounds for missing values.}
  \item{tab.u}{A numeric matrix of upper bounds for missing values.}
  \item{res}{An output list resulting from the function \code{\link{estim.mix}}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric matrix of estimated probabilities to be MCAR for missing values in the confidence intervals defined thanks to \code{tab.l} and \code{tab.u}.
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

#Imputation of missing values with the slsa algorithm
dat.slsa=impute.slsa(tab=res.sim$dat.obs,conditions=res.sim$condition,repbio=res.sim$repbio);

#Estimation of the mixture model
res=estim.mix(tab=res.sim$dat.obs, tab.imp=dat.slsa, conditions=res.sim$condition);

#Computing probabilities to be MCAR
born=estim.bound(tab=res.sim$dat.obs,conditions=res.sim$condition);
proba=prob.mcar.tab(born$tab.lower,born$tab.upper,res);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
