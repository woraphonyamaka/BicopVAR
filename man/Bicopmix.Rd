\name{Bicopmix}
\alias{Bicopmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Bivariate Mix-copula based Vector autoregressive model
}
\description{
This function introduces a mixed copula-based vector autoregressive (VAR) model for investigating
the relationship between random variables. The one-step maximum likelihood estimation is used to obtain point estimates of the autoregressive parameters and mixed copula parameters. More specifically, we combine the likelihoods of the marginal and mixed copula to construct the full likelihood function.
}
\usage{
Bicopmix(Y, p,dist1,dist2,copula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y}{ Vector of dependent variable}
\item{family}{ Copula functions eg. normal=normal,Student-t=studentt,clayton=clayton,gumbel=gumbel,
                frank=frank,joe=joe,Normal-Studenet-t=NS,Normal-clayton=NC,clayton-frank=CF, ..., etc.}
\item{dist1}{  Margninal distribution of first endogenous variable}
\item{dist2}{  Margninal distribution of second endogenous variable}
\item{p}{ The numbe of lag term}
}
\details{
The one step estimation of full conditional likelihood function of Bivariate Mix-copula based Vector autoregressive model is used. Several single and mixed copula functions are introduced.
}
\value{
\item{result }{The result contain the estimated parameters, standard errors, t-stat, and p-value}
\item{AIC }{Akaiki Information Criteria}
\item{BIC }{Bayesian Information Criteria}
\item{Loglikelihood }{Maximum Log-likelihood function}
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

# simulation
n=200
mC <- mixCopula(list(gumbelCopula(2.5, dim=2),
claytonCopula(1.8, dim=2)),c(0.2,0.8))

set.seed(17)
uM <- rCopula(n, mC)
e1=qnorm(uM[,1],0,1.5)
e2=qnorm(uM[,2],0,1)
cor(cbind(e1,e2))
e.var = cbind(e1,e2)
e.var = t(e.var)
pi1 = matrix(c(0.7,0.2,0.2,0.7),2,2)
mu.vec = c(1,5)
c.vec = as.vector((diag(2)-pi1)%*%mu.vec)
y.var = matrix(0,2,n)
y.var[,1] = mu.vec
for (i in 2:n) {
     y.var[,i] = c.vec+pi1%*%y.var[,i-1]+e.var[,i]
}
Y= t(y.var)

## Estimation here
p=1
m=ncol(Y)
dist1="norm"
dist2="norm"
copula="GJ"
fit=Bicopmix(Y,p,dist1,dist2,copula)
fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mixed Copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ VAR }% __ONLY ONE__ keyword per line
