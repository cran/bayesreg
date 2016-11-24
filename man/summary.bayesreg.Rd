\name{summary.bayesreg}
\alias{summary.bayesreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize Bayesian penalised regression (bayesreg) fits
}
\description{
\code{summary} method for class \code{"bayesreg"}.
}
\usage{
\method{summary}{bayesreg}(object, sortrank = FALSE, displayor = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{"bayesreg"} created as a result of a call to \code{\link{bayesreg}}.
}
  \item{sortrank}{logical; if \code{TRUE}, the variables in the summary will be sorted by their importance as determined by their rank estimated by the Bayesian feature ranking algorithm.}
  
  \item{displayor}{logical; if \code{TRUE}, the variables will be summarised in terms of their cross-sectional odds-ratios rather than their regression coefficients.}
  
  \item{...}{further arguments passed to or from other methods.}
}
\details{
The \code{summary} method computes a number of summary statistics and displays these for each variable in a table, along with suitable header information.

For continuous target variables, the header information includes an estimate of the residual squared error, the $R^2$ statistic and the Deviance information criterion (DIC) statistic. For logistic regression models, the header information includes the negative log-likelihood at the posterior mean of the regression coefficients, the pseudo $R^2$ score and the DIC statistic.

The main table summarises properties of the coefficients for each of the variables. The first column is the variable name. The second and third columns are either the mean and standard error of the coefficients, or the median and standard error of the cross-sectional odds-ratios if \code{displayor=TRUE}. 

The fourth and fifth columns are the 95\% credible intervals of the coefficients (odds-ratios). The sixth column displays the posterior $t$-statistic, calculated as the ratio of the posterior mean on the posterior standard deviation for the coefficient. The seventh column is the importance rank assigned to the variable by the Bayesian feature ranking algorithm. 

In betwen the seventh and eighth columns are up to two asterisks indicating significance; a variable scores a first asterisk if the 75\% credible interval does not include zero, and scores a second if the 95\% credible interval does not include zero. The final column gives an estimate of the effective sample size for the variable, ranging from 0 to 100, which indicates how much autocorrelation is present in the sample chain for this variable (smaller values indicating high levels of autocorrelation).
}
\value{
Returns an object with the following fields:

\item{logl}{The log-likelihood of the model at the posterior mean estimates of the regression coefficients.}
\item{dic}{The deviance information criterion score of the model.}
\item{rsqr}{For non-binary data, the R^2 statistic.}
\item{rootmse}{For non-binary data, the estimated root-mean-squared residual error.}
\item{p_rsqr}{For binary data, the pseudo-R^2 statistic.}
\item{muCoef}{The posterior means of the regression coefficients.}
\item{seCoef}{The posterior standard deviations of the regression coefficients.}
\item{CICoef}{The posterior 95\% credible interval for the regression coefficients.}
\item{medOR}{For binary data, the posterior median of the cross-sectional odds-ratios.}
\item{seOR}{For binary data, the posterior standard deviation of the cross-sectional odds-ratios.}
\item{CIOR}{For binary data, the posterior 95\% credible interval for the cross-sectional odds-ratios.}
\item{tStat}{The posterior t-statistic for the coefficients.}
\item{nStars}{The significance level for the variable (see above).}
\item{rank}{The variable importance rank as estimated by the Bayesian feature ranking algorithm (see above).}
\item{ESS}{The effective sample size for the variable.}

}
\author{
Enes Makalic and Daniel F. Schmidt

Maintainer: Daniel F. Schmidt \email{dschmidt@unimelb.edu.au}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The model fitting function \code{\link{bayesreg}} and prediction function \code{\link{predict.bayesreg}}
}
\examples{

X = matrix(rnorm(100*20),100,20)
b = matrix(0,20,1)
b[1:9] = c(0,0,0,0,5,4,3,2,1)
y = X \%*\% b + rnorm(100, 0, 1)
df <- data.frame(X,y)

rv_hs <- bayesreg(y~.,df,prior="hs")       # Horseshoe regression

# Summarise without sorting by variable rank
rv_hs.s <- summary(rv_hs)

# Summarise sorting by variable rank
rv_hs.s <- summary(rv_hs, sortrank = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
