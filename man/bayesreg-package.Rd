% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesreg-package.R
\docType{package}
\name{bayesreg-package}
\alias{bayesreg-package}
\title{Getting started with the Bayesreg package}
\description{
This is a comprehensive, user-friendly package implementing the state-of-the-art in 
Bayesian linear regression, Bayesian count regression and Bayesian logistic regression. Features of the toolbox include:
\itemize{
 \item Supports Gaussian, Laplace, Student-t, Poisson, geometric and logistic binary data models.
 \item Efficient and numerically stable implementations of Bayesian ridge, Bayesian lasso, horseshoe and horseshoe+ regression.
 \item Provides variable ranking and importance, credible intervals and diagnostics such as the widely applicable information criterion.
 \item Factor variables are automatically grouped together and additional shrinkage is applied to the set of indicator variables to which they expand.
 \item Prediction tools for generating credible intervals and Bayesian averaging of predictions.
 \item Support for multiple cores
}
The lasso, horseshoe and horseshoe+ priors are recommended for data sets where the number of 
predictors is greater than the sample size. The Laplace, Student-t and logistic models are based on scale-mixture representations;
logistic regression utilises the Polya-gamma sampler implemented in the \code{pgdraw} package. The Poisson and
geometric distributions are implemented using a fast gradient-assisted Metropolis-Hastings algorithm.
}
\details{
Count (non-negative integer) regression is now supported through implementation of Poisson and geometric regression models.
To support analysis of data with outliers, we provide two heavy-tailed error models in our 
implementation of Bayesian linear regression: Laplace and Student-t distribution errors. 
The widely applicable information criterion (WAIC) is routinely calculated and displayed
to assist users in selecting an appropriate prior distribution for their particular problem, i.e., choice of regularisation or data model.
Most features are straightforward to use. The package will make use of multiple CPU cores by default, if available. 
This feature may be disabled if necessary for problems with very large predictor matrices.

Further information on the particular algorithms/methods implemented in this package provided
by the literature referenced below.

Version history:
\itemize{
\item Version 1.1: Initial release.
\item Version 1.2: Added Poisson and geometric regression; user specifiable credible interval levels for \code{summary()} and \code{predict()}; \code{summary()} column "ESS" now reports effective sample size rather than percentage-effective sample size.
\item Version 1.3: Added support for multiple cores; if n.cores = Inf (the default) Bayesreg will divide the requested number of samples by the number of available cores (total cores minus one) and run these as seperate sampling chains, in parallel, and then recombine after sampling. It is possible to explicitly control the number of cores being used via the n.cores option if desired; n.cores=1 will disable parallelization.
}
}
\note{
To cite this package please reference: 

Makalic, E. & Schmidt, D. F.
High-Dimensional Bayesian Regularised Regression with the BayesReg Package
arXiv:1611.06649 [stat.CO], 2016 \url{http://arxiv.org/pdf/1611.06649}

A MATLAB-compatible implementation of this package can be obtained from:

\url{https://au.mathworks.com/matlabcentral/fileexchange/60823-flexible-bayesian-penalized-regression-modelling}
}
\examples{
\dontrun{
# -----------------------------------------------------------------
# By default Bayesreg now utilizes multiple cores if you have them
# available. If you do not want to use multiple cores you can set 
# n.cores=1 when calling Bayesreg.
#
# In most realistic/practical settings, i.e., when a 
# even a moderate number of samples is being requested, parallelization 
# will usually result in substantial speedups, and can dramatically reduce
# run-time for large numbers of samples.
#
# However, if your design matrix and desired number of samples 
# are small (i.e, small nrow(X), ncol(X) and n.samples) 
# then sometimes no parallelization, or parallelization with a small
# number of cores, can be quicker due to the overhead of setting up 
# multiple threads. This is likely only relevant if you are calling Bayesreg
# hundreds/thousands of times with *small* numbers of samples being requested
# for each call, e.g., if you are doing some sort of one-at-a-time testing
# across many predictors such as a genome-wide association test, or similar.
# In this case you may be better off parallelizing the calls to 
# Bayesreg across the tests and disabling parallelization when calling 
# Bayesreg (i.e., n.cores=1).
#
# -----------------------------------------------------------------
# Example 1: Gaussian regression
#
X = matrix(rnorm(100*20),100,20)
b = matrix(0,20,1)
b[1:5] = c(5,4,3,2,1)
y = X \%*\% b + rnorm(100, 0, 1)

df <- data.frame(X,y)
rv.lm <- lm(y~.,df)                        # Regular least-squares
summary(rv.lm)

# Horseshoe regression -- here we show how to explicitly control the maximum
# number of cores being used for sampling to 4
rv.hs <- bayesreg(y~., df, prior="hs", n.cores=4)       
rv.hs$n.cores  # actual number of cores used (will be <= 4, depending on machine)
rv.hs.s <- summary(rv.hs)

# Expected squared prediction error for least-squares
coef_ls = coef(rv.lm)
as.numeric(sum( (as.matrix(coef_ls[-1]) - b)^2 ) + coef_ls[1]^2)

# Expected squared prediction error for horseshoe
as.numeric(sum( (rv.hs$mu.beta - b)^2 ) + rv.hs$mu.beta0^2)


# -----------------------------------------------------------------
# Example 2: Gaussian v Student-t robust regression
X = 1:10;
y = c(-0.6867, 1.7258, 1.9117, 6.1832, 5.3636, 7.1139, 9.5668, 10.0593, 11.4044, 6.1677);
df = data.frame(X,y)

# Gaussian ridge
rv.G <- bayesreg(y~., df, model = "gaussian", prior = "ridge", n.samples = 1e5)

# Student-t ridge
rv.t <- bayesreg(y~., df, model = "t", prior = "ridge", t.dof = 5, n.samples = 1e5)

# Plot the different estimates with credible intervals
plot(df$X, df$y, xlab="x", ylab="y")

yhat_G <- predict(rv.G, df, bayes.avg=TRUE)
lines(df$X, yhat_G[,1], col="blue", lwd=2.5)
lines(df$X, yhat_G[,3], col="blue", lwd=1, lty="dashed")
lines(df$X, yhat_G[,4], col="blue", lwd=1, lty="dashed")

yhat_t <- predict(rv.t, df, bayes.avg=TRUE)
lines(df$X, yhat_t[,1], col="darkred", lwd=2.5)
lines(df$X, yhat_t[,3], col="darkred", lwd=1, lty="dashed")
lines(df$X, yhat_t[,4], col="darkred", lwd=1, lty="dashed")

legend(1,11,c("Gaussian","Student-t (dof=5)"),lty=c(1,1),col=c("blue","darkred"),
       lwd=c(2.5,2.5), cex=0.7)


# -----------------------------------------------------------------
# Example 3: Poisson/geometric regression example

X  = matrix(rnorm(100*5),100,5)
b  = c(0.5,-1,0,0,1)
nu = X\%*\%b + 1
y  = rpois(lambda=exp(nu),n=length(nu))

df <- data.frame(X,y)

# Fit a Poisson regression
rv.pois=bayesreg(y~., data=df, model="poisson",prior="hs", burnin=1e4, n.samples=5e4)
summary(rv.pois)

# Fit a geometric regression
rv.geo=bayesreg(y~., data=df, model="geometric",prior="hs", burnin=1e4, n.samples=5e4)
summary(rv.geo)

# Compare the two models in terms of their WAIC scores
cat(sprintf("Poisson regression WAIC=\%g vs geometric regression WAIC=\%g", 
            rv.pois$waic, rv.geo$waic))
# Poisson is clearly preferred to geometric, which is good as data is generated from a Poisson!
 
 
# -----------------------------------------------------------------
# Example 4: Logistic regression on spambase
data(spambase)
  
# bayesreg expects binary targets to be factors
spambase$is.spam <- factor(spambase$is.spam)

# First take a subset of the data (1/10th) for training, reserve the rest for testing
spambase.tr  = spambase[seq(1,nrow(spambase),10),]
spambase.tst = spambase[-seq(1,nrow(spambase),10),]
  
# Fit a model using logistic horseshoe for 2,000 samples
# In practice, > 10,000 samples would be a more realistic amount to draw
rv <- bayesreg(is.spam ~ ., spambase.tr, model = "logistic", prior = "horseshoe", n.samples = 2e3)
  
# Summarise, sorting variables by their ranking importance
rv.s <- summary(rv,sort.rank=TRUE)
  
# Make predictions about testing data -- get class predictions and class probabilities
y_pred <- predict(rv, spambase.tst, type='class')
  
# Check how well did our predictions did by generating confusion matrix
table(y_pred, spambase.tst$is.spam)
  
# Calculate logarithmic loss on test data
y_prob <- predict(rv, spambase.tst, type='prob')
cat('Neg Log-Like for no Bayes average, posterior mean estimates: ', sum(-log(y_prob[,1])), '\n')
y_prob <- predict(rv, spambase.tst, type='prob', sum.stat="median")
cat('Neg Log-Like for no Bayes average, posterior median estimates: ', sum(-log(y_prob[,1])), '\n')
y_prob <- predict(rv, spambase.tst, type='prob', bayes.avg=TRUE)
cat('Neg Log-Like for Bayes average: ', sum(-log(y_prob[,1])), '\n')
}

}
\references{
Bhadra, A.; Datta, J.; Polson, N. G. & Willard, B. 
   The Horseshoe+ Estimator of Ultra-Sparse Signals 
   Bayesian Analysis, 2016

   Bhattacharya, A.; Chakraborty, A. & Mallick, B. K. 
   Fast sampling with Gaussian scale-mixture priors in high-dimensional regression 
   arXiv:1506.04778, 2016

   Carvalho, C. M.; Polson, N. G. & Scott, J. G. 
   The horseshoe estimator for sparse signals 
   Biometrika, Vol. 97, pp. 465-480, 2010

   Makalic, E. & Schmidt, D. F. 
   A Simple Sampler for the Horseshoe Estimator 
   IEEE Signal Processing Letters, Vol. 23, pp. 179-182, 2016

   Park, T. & Casella, G. 
   The Bayesian Lasso 
   Journal of the American Statistical Association, Vol. 103, pp. 681-686, 2008

   Polson, N. G.; Scott, J. G. & Windle, J. 
   Bayesian inference for logistic models using Polya-Gamma latent variables 
   Journal of the American Statistical Association, Vol. 108, pp. 1339-1349, 2013

   Rue, H. 
   Fast sampling of Gaussian Markov random fields 
   Journal of the Royal Statistical Society (Series B), Vol. 63, pp. 325-338, 2001

   Xu, Z., Schmidt, D.F., Makalic, E., Qian, G. & Hopper, J.L.
   Bayesian Grouped Horseshoe Regression with Application to Additive Models
   AI 2016: Advances in Artificial Intelligence, pp. 229-240, 2016

   Schmidt, D.F. & Makalic, E.
   Bayesian Generalized Horseshoe Estimation of Generalized Linear Models
   ECML PKDD 2019: Machine Learning and Knowledge Discovery in Databases. pp 598-613, 2019
   
   Stan Development Team, Stan Reference Manual (Version 2.26), Section 15.4, "Effective Sample Size",
   \url{https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html}
}
\seealso{
\code{\link{bayesreg}}
}
\author{
Daniel Schmidt \email{daniel.schmidt@monash.edu} 

Department of Data Science and AI, Monash University, Australia

Enes Makalic \email{enes.makalic@monash.edu}

Department of Data Science and AI, Monash University, Australia
}
