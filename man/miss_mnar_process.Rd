\name{miss.mnar.process}
\alias{miss.mnar.process}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the MNAR mechanism in a sample.
}
\description{
This function allows estimating the MNAR data mechanism, i.e. the probability to be MNAR given that the value is missing in function of the intensity level, from an estimation of a mixture model of MNAR and MCAR values (see \code{\link{estim.mix}} function).
}
\usage{
miss.mnar.process(abs,pi_mcar,F_mnar,F_na)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{abs}{The interval on which is estimated the MNAR data mechanism.}
    \item{pi_mcar}{An estimation of the proportion of MCAR values.}
  \item{F_mnar}{An estimation of the cumulative distribution function of the MNAR values on the interval \code{abs}.}
\item{F_na}{An estimation of the cumulative distribution function of the missing values on the interval \code{abs}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list composed of:
\item{abs}{The interval on which is estimated the MNAR data mechanism.}
\item{p}{The estimated probability to be MNAR given that the value is missing in function of the intensity level.}
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

#Estimating the MNAR mechanism in the first replicate
mnp=miss.mnar.process(res$abs.mod,res$pi.mcar[1],res$F.mnar[,1],res$F.na[,1])
plot(mnp$abs,mnp$p,ty="l",xlab="Intensity values",ylab="Estimated probability to be MNAR")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
