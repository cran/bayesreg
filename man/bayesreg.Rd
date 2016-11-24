\name{bayesreg}
\alias{bayesreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Bayesian Regression Models with Continuous Shrinkage Priors
}
\description{
Fit a linear or logistic regression model using Bayesian continuous shrinkage prior distributions. Handles ridge, lasso, horseshoe and horseshoe+ regression with logistic, Gaussian, Laplace or Student-t distributed targets.
}
\usage{
bayesreg(formula, data, model = "normal", prior = "ridge",  
     nsamples = 1000, burnin = 1000, thin = 5, tdof = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "\code{\link{formula}}": a symbolic description of the model to be fitted using the standard R formula notation.
%%     ~~Describe \code{formula} here~~
}
  \item{data}{an data frame containing the variables in the model.
%%     ~~Describe \code{data} here~~
}
  \item{model}{the distribution of the target (y) variable. Continuous or numeric variables can be distributed as per a Gaussian distribution (\code{model="gaussian"} or \code{model="normal"}), Laplace distribution (\code{model = "laplace"} or \code{model = "l1"}) or Student-t distribution (\code{"model" = "studentt"} or \code{"model" = "t"}). For binary targets (factors with two levels) either \code{model="logistic"} or \code{"model"="binomial"} should be used.
%%     ~~Describe \code{model} here~~
}
  \item{prior}{The continuous shrinkage prior distribution over the regression coefficients that is to be used. Options include ridge regression (\code{prior="rr"} or \code{prior="ridge"}), lasso regression (\code{prior="lasso"}), horseshoe regression (\code{prior="hs"} or \code{prior="horseshoe"}) and horseshoe+ regression (\code{prior="hs+"} or \code{prior="horseshoe+"}).
%%     ~~Describe \code{prior} here~~
}
  \item{nsamples}{number of posterior samples to generate.}
  \item{burnin}{number of burn-in samples.}
  \item{thin}{level of thinning.}
  \item{tdof}{degrees of freedom for the Student-t distribution.}
}
\details{
Draws a series of samples from the posterior distribution of a linear (Gaussian, Laplace or Student-t) or logistic regression model with specified continuous shrinkage prior distribution (ridge regression, lasso, horseshoe and horseshoe+) using Gibbs sampling. The intercept parameter is always included, and is never penalised.

While only \code{nsamples} are returned, the total number of samples generated is equal to \code{burnin}+\code{nsamples}*\code{thin}. To generate the samples of the regression coefficients, the code will use either Rue's algorithm (when the number of samples is twice the number of covariates) or the algorithm of Bhattacharya et al. as appropriate.
}
\value{
An object with S3 class \code{"bayesreg"}. 

\item{beta}{Posterior samples the regression model coefficients}
\item{beta0}{Posterior samples of the intercept parameter}
\item{sigma2}{Posterior samples of the square of the scale parameter; for Gaussian distributed targets this is equal to the variance. For binary distributed targets this is empty.}
\item{muBeta}{The mean of the posterior samples for the regression coefficients.}
\item{muBeta0}{The mean of the posterior samples for the intercept parameter.}
\item{muSigma2}{The mean of the posterior samples for squared scale parameter.}

\item{tau2}{Posterior samples of the global shrinkage parameter.}

\item{tstat}{Posterior t-statistics for each regression coefficient.}

\item{varranks}{Ranking of the covariates by their importance, with "1" denoting the most important covariate.}

\item{logl}{The log-likelihood at the posterior means of the model parameters}

\item{dic}{The Deviance Information Criterion (DIC) score for the model}

The returned object also stores the parameters/options used to run \code{bayesreg}:

\item{formula}{The object of type "\code{\link{formula}}" describing the fitted model.}
\item{model}{The distribution of the target (y) variable.}
\item{prior}{The shrinkage prior used to fit the model.}
\item{nsamples}{The number of samples generated from the posterior distribution.}
\item{burnin}{The number of burnin samples that were generated.}
\item{thin}{The level of thinning.}
\item{n}{The sample size of the data used to fit the model.}
\item{p}{The number of covariates in the fitted model.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Enes Makalic and Daniel F. Schmidt

Maintainer: Daniel F. Schmidt \email{dschmidt@unimelb.edu.au}
}
\references{
      Makalic, E. & Schmidt, D. F.
      High-Dimensional Bayesian Regularised Regression with the BayesReg Package
      arXiv:1611.06649 [stat.CO], 2016 \url{https://arxiv.org/pdf/1611.06649.pdf}
 
      Park, T. & Casella, G. 
      The Bayesian Lasso 
      Journal of the American Statistical Association, Vol. 103, pp. 681-686, 2008
 
      Carvalho, C. M.; Polson, N. G. & Scott, J. G. 
      The horseshoe estimator for sparse signals 
      Biometrika, Vol. 97, 465-480, 2010
 
      Makalic, E. & Schmidt, D. F. 
      A Simple Sampler for the Horseshoe Estimator 
      IEEE Signal Processing Letters, Vol. 23, pp. 179-182, 2016 \url{https://arxiv.org/pdf/1508.03884v4.pdf}
 
      Bhadra, A.; Datta, J.; Polson, N. G. & Willard, B. 
      The Horseshoe+ Estimator of Ultra-Sparse Signals 
      Bayesian Analysis, 2016
 
      Polson, N. G.; Scott, J. G. & Windle, J. 
      Bayesian inference for logistic models using Polya-Gamma latent variables 
      Journal of the American Statistical Association, Vol. 108, 1339-1349, 2013
 
      Rue, H. 
      Fast sampling of Gaussian Markov random fields 
      Journal of the Royal Statistical Society (Series B), Vol. 63, 325-338, 2001
 
      Bhattacharya, A.; Chakraborty, A. & Mallick, B. K. 
      Fast sampling with Gaussian scale-mixture priors in high-dimensional regression 
      arXiv:1506.04778, 2016
}
\note{
    To cite this toolbox please reference: 
    
      Makalic, E. & Schmidt, D. F.
      High-Dimensional Bayesian Regularised Regression with the BayesReg Package
      arXiv:1611.06649 [stat.CO], 2016 \url{https://arxiv.org/pdf/1611.06649.pdf}
      
      A MATLAB implementation of the bayesreg function is also available from:
      
      \url{https://au.mathworks.com/matlabcentral/fileexchange/60335-bayesian-regularised-linear-and-logistic-regression}
}
      
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predict.bayesreg}} and \code{\link{summary.bayesreg}}
}
\examples{

# -----------------------------------------------------------------
# Example 1: Gaussian regression
X = matrix(rnorm(100*20),100,20)
b = matrix(0,20,1)
b[1:5] = c(5,4,3,2,1)
y = X \%*\% b + rnorm(100, 0, 1)

df <- data.frame(X,y)
rv_lm <- lm(y~.,df)                        # Regular least-squares
summary(rv_lm)

rv_hs <- bayesreg(y~.,df,prior="hs")       # Horseshoe regression
rv_hs.s <- summary(rv_hs)

# Expected squared prediction error for least-squares
coef_ls = coef(rv_lm)
as.numeric(sum( (as.matrix(coef_ls[-1]) - b)^2 ) + coef_ls[1]^2)

# Expected squared prediction error for horseshoe
as.numeric(sum( (rv_hs$muBeta - b)^2 ) + rv_hs$muBeta0^2)


# -----------------------------------------------------------------
# Example 2: Gaussian v Student-t robust regression
X = 1:10;
y = c(-0.6867, 1.7258, 1.9117, 6.1832, 5.3636, 7.1139, 9.5668, 10.0593, 11.4044, 6.1677);
df = data.frame(X,y)

# Gaussian ridge
rv_G <- bayesreg(y~., df, model = "gaussian", prior = "ridge", nsamples = 1e3)

# Student-t ridge
rv_t <- bayesreg(y~., df, model = "t", prior = "ridge", tdof = 5, nsamples = 1e3)

# Plot the different estimates with credible intervals
plot(df$X, df$y, xlab="x", ylab="y")

yhat_G <- predict(rv_G, df, bayesavg=TRUE)
lines(df$X, yhat_G[,1], col="blue", lwd=2.5)
lines(df$X, yhat_G[,3], col="blue", lwd=1, lty="dashed")
lines(df$X, yhat_G[,4], col="blue", lwd=1, lty="dashed")

yhat_t <- predict(rv_t, df, bayesavg=TRUE)
lines(df$X, yhat_t[,1], col="darkred", lwd=2.5)
lines(df$X, yhat_t[,3], col="darkred", lwd=1, lty="dashed")
lines(df$X, yhat_t[,4], col="darkred", lwd=1, lty="dashed")

legend(1,11,c("Gaussian","Student-t (dof=5)"),lty=c(1,1),col=c("blue","darkred"),
       lwd=c(2.5,2.5), cex=0.7)

\dontrun{
# -----------------------------------------------------------------
# Example 3: Logistic regression on spambase
data(spambase)

# bayesreg expects binary targets to be factors
spambase$is.spam <- factor(spambase$is.spam)

# First take a subset of the data (1/10th) for training, reserve the rest for testing
spambase.tr  = spambase[seq(1,nrow(spambase),10),]
spambase.tst = spambase[-seq(1,nrow(spambase),10),]

# Fit a model using logistic horseshoe for 2,000 samples
rv <- bayesreg(is.spam ~ ., spambase.tr, model = "logistic", prior = "horseshoe", nsamples = 2e3)

# Summarise, sorting variables by their ranking importance
rv.s <- summary(rv,sortrank=TRUE)

# Make predictions about testing data -- get class predictions and class probabilities
y_pred <- predict(rv, spambase.tst, type='class')

# Check how well did our predictions did by generating confusion matrix
table(y_pred, spambase.tst$is.spam)

# Calculate logarithmic loss on test data
y_prob <- predict(rv, spambase.tst, type='prob')
cat('Neg Log-Like for no Bayes average, posterior mean estimates: ', sum(-log(y_prob[,1])))
y_prob <- predict(rv, spambase.tst, type='prob', sum.stat="median")
cat('Neg Log-Like for no Bayes average, posterior median estimates: ', sum(-log(y_prob[,1])))
y_prob <- predict(rv, spambase.tst, type='prob', bayesavg=TRUE)
cat('Neg Log-Like for Bayes average: ', sum(-log(y_prob[,1])))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
