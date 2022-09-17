\name{IRF}
\alias{IRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impulse response function of Bivariate Mix-copula based Vector autoregressive model
}
\description{
An impulse-response function describes the evolution of the variable of interest along a
specified time horizon after a shock in a given moment
}
\usage{
IRF=function(Y,A,Sigma,h.step)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y}{ Vector of dependent variable}
\item{A}{ 2x2xp Autoregressive coefficeint matrix}
\item{Sigma}{  2x2 Variance covariance matrix }
\item{h.step}{ Integer specifying the steps}
}
\details{
Computes the impulse response coefficients of a Bivariate Mixed Copula based VAR(p) for h.steps aheads.
}
\value{
\item{IRFplot}{ Impulse response coefficients plots}
}
\references{
Yamaka, W., & Thongkairat, S. (2020). A Mixed Copula-Based Vector Autoregressive Model for Econometric Analysis. International Journal of Uncertainty, Fuzziness and Knowledge-Based Systems, 28(Supp01), 113-121.

Yamaka, W., Gupta, R., Thongkairat, S., & Maneejuk, P. (2021). Structural and Predictive Analyses with a Mixed Copula-Based Vector Autoregression Model, Journal of Forecasting, https://doi.org/10.1002/for.2902

}
\author{
Woraphon Yamaka, Paravee Maneejuk, and Sukrit Thongkairat
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Required packages
library(mvtnorm)
library(vars)
library(copula)
library(VineCopula)

#=== Impulse response function
1 = Gaussian copula
2 = Student t copula (t-copula)
3 = Clayton copula
4 = Gumbel copula
5 = Frank copula
6 = Joe copula
K=length(fit$result[,1])
sigma1=fit$result[K-2,1]
sigma2=fit$result[K-1,1]
tau1 <- BiCopPar2Tau(family = 5, par = sigma1)
tau2 <- BiCopPar2Tau(family = 6, par = sigma2)
weight=fit$result[K]
Etau=(weight*tau2)+((1-weight)*tau1)
covsigma=Etau*sqrt(fit$result[K-4,1])*sqrt(fit$result[K-3,1])
Phi=fit$result[c(2,3,5,6),1]
A=matrix(Phi,2,2)
Sigma=matrix(c(sigma1,covsigma,covsigma,sigma2),2,2)
IRF(Y,A,Sigma,h.step=10)

\keyword{ Impulse response  }% use one of  RShowDoc("KEYWORDS")
\keyword{ VAR }% __ONLY ONE__ keyword per line
