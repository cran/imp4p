\name{delete.na.rows}
\alias{delete.na.rows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Deleting the rows of an input matrix with only missing values in a condition after using \code{\link{sim.data}}.
}
\description{
This function allows to delete rows without any observed values in one of all the conditions. It is built to be used after the function \code{\link{sim.data}}.
}
\usage{
delete.na.rows(tab,tab.c,conditions,list.MCAR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to the intensities measured in an experimental sample, and each row to the ones of an identified peptide.}
    \item{tab.c}{A matrix where the missing values of \code{tab} have been imputed under the assumption that they are all MCAR. For instance, such a matrix can be obtained by using the function \code{\link{impute.slsa}} of this package.}
  \item{conditions}{A vector of factors indicating the biological condition to which each sample belongs.}
  \item{list.MCAR}{A matrix containing the indexes of the rows in \code{tab} which contain MCAR values in each column.}
}
\details{
All the rows of an input matrix are deleted if they have only missing values in one of all the conditions. This function must be used after the function \code{\link{sim.data}} to keep only rows with at least one observed value in each condition.
}
\value{
A list composed of:
  \item{tab.mod}{The matrix \code{tab} with deleted rows.}
  \item{tab.comp}{The matrix \code{tab.c} with deleted rows.}
  \item{list.delete}{A vector containing the indexes of the rows in \code{tab} which have been deleted.}
  \item{pi.na}{Proportion of missing values in each column of \code{tab.mod}.}
  \item{pi.mcar}{Proportion of MCAR values in each column of \code{tab.mod}.}
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
\code{
\link{sim.data}
}
}

\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600,pi.mcar=0.2,para=10,nb.cond=2,nb.repbio=3,
nb.sample=5,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2);

#Deleting rows without any observed value in a condition
result=delete.na.rows(tab=res.sim$dat.obs, tab.c=res.sim$dat.comp, conditions=res.sim$conditions,
list.MCAR=res.sim$list.MCAR);

#Determining the lists of MNAR and MCAR values in the first column of the dataset
#after deletion of rows
listmcar=which((1:length(res.sim$dat.obs[,1])\%in\%res.sim$list.MCAR[,1])[-result$list.delete]);

listmnar=which(is.na(result$tab.mod[,1]))[!which(is.na(result$tab.mod[,1]))\%in\%listmcar];
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
