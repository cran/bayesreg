\name{predict.bayesreg}
\alias{predict.bayesreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Bayesian penalised regression (bayesreg) fits
}
\description{
Predict values based on Bayesian penalised regression fits
}
\usage{
\method{predict}{bayesreg}(object, df, type = "linpred", bayesavg = FALSE, sum.stat = "mean", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"bayesreg"} created as a result of a call to \code{\link{bayesreg}}.}
%%     ~~Describe \code{object} here~~

  \item{df}{A data frame providing the variables from which to produce predictions.
%%     ~~Describe \code{\dots} here~~
}
  \item{type}{The type of predictions to produce; if \code{type="linpred"} it will return the linear predictor for both binary and continuous data. For binary data, if \code{type="prob"} it will return predictive probability estimates, and if \code{type="class"} and the data is binary, it will return the best guess at the class of the target variable.}
  
  \item{bayesavg}{Whether to produce predictions using Bayesian averaging.}
  
  \item{sum.stat}{The type of summary statistic to use; either \code{sum.stat="mean"} or \code{sum.stat="median"}.}
  
  \item{...}{further arguments passed to or from other methods.}
}
\details{
\code{predict.bayesreg} produces predicted values using variables from the specified data frame. The type of predictions produced depend on the value of the parameter \code{type}.

If \code{type="linpred"}, the predictions that are returned will be the value of the linear predictor formed from the model coefficients and the provided data. 

If \code{type="prob"}, the predictions will be probabilities. If the specified data frame includes a column with the same name as the target variable on which the model was created, the predictions will then be the probability density values for these target values. For binary data, if the specified data frame does not include a column with the same name as the target variable, the predictions will be the probability of the target being equal to the second level of the factor variable. 

If \code{type="class"} and the target variable is binary, the predictions will be the most likely class.

If \code{bayesavg} is \code{FALSE} the predictions will be produced by using a summary of the posterior samples of the coefficients and scale parameters as estimates for the model. If \code{bayesavg} is \code{TRUE}, the predictions will be produced by posterior averaging over the posterior samples of the coefficients and scale parameters, allowing the uncertainty in the estimation process to be explicitly taken into account in the prediction process. 

If \code{sum.stat="mean"} and \code{bayesavg} is \code{FALSE}, the mean of the posterior samples will be used as point estimates for making predictions. Likewise, if \code{sum.stat="median"} and \code{bayesavg} is \code{FALSE}, the co-ordinate wise posterior medians will be used as estimates for making predictions. If \code{bayesavg} is \code{TRUE} and \code{type!="prob"}, the posterior mean (median) of the predictions from each of the posterior samples will be used as predictions. The value of \code{sum.stat} has no effect if \code{type="prob"}.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{predict.bayesreg} produces a vector or matrix of predictions of the specified type. If \code{bayesavg} is \code{FALSE} a matrix with a single column \code{pred} is returned, containing the predictions.

If \code{bayesavg} is \code{TRUE}, three additional columns are returned: \code{se(pred)}, which contains standard errors for the predictions, and \code{CI 5\%} and \code{CI 95\%} which contain 95\% credible intervals for the predictions.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
The model fitting function \code{\link{bayesreg}} and summary function \code{\link{summary.bayesreg}}
}
\examples{

# -----------------------------------------------------------------
# Example 1: Fitting linear models to data and generating credible intervals
X = 1:10;
y = c(-0.6867, 1.7258, 1.9117, 6.1832, 5.3636, 7.1139, 9.5668, 10.0593, 11.4044, 6.1677);
df = data.frame(X,y)

# Gaussian ridge
rv_L <- bayesreg(y~., df, model = "laplace", prior = "ridge", nsamples = 1e3)

# Plot the different estimates with credible intervals
plot(df$X, df$y, xlab="x", ylab="y")

yhat <- predict(rv_L, df, bayesavg=TRUE)
lines(df$X, yhat[,1], col="blue", lwd=2.5)
lines(df$X, yhat[,3], col="blue", lwd=1, lty="dashed")
lines(df$X, yhat[,4], col="blue", lwd=1, lty="dashed")
yhat <- predict(rv_L, df, bayesavg=TRUE, sum.stat = "median")
lines(df$X, yhat[,1], col="red", lwd=2.5)

legend(1,11,c("Posterior Mean (Bayes Average)","Posterior Median (Bayes Average)"),
       lty=c(1,1),col=c("blue","red"),lwd=c(2.5,2.5), cex=0.7)


# -----------------------------------------------------------------
# Example 2: Predictive density for continuous data
X = 1:10;
y = c(-0.6867, 1.7258, 1.9117, 6.1832, 5.3636, 7.1139, 9.5668, 10.0593, 11.4044, 6.1677);
df = data.frame(X,y)

# Gaussian ridge
rv_G <- bayesreg(y~., df, model = "gaussian", prior = "ridge", nsamples = 1e3)

# Produce predictive density for X=2
df.tst = data.frame(y=seq(-7,12,0.01),X=2)
prob_noavg_mean <- predict(rv_G, df.tst, bayesavg=FALSE, type="prob", sum.stat = "mean")
prob_noavg_med  <- predict(rv_G, df.tst, bayesavg=FALSE, type="prob", sum.stat = "median")
prob_avg        <- predict(rv_G, df.tst, bayesavg=TRUE, type="prob")

# Plot the density
plot(NULL, xlim=c(-7,12), ylim=c(0,0.14), xlab="y", ylab="p(y)")
lines(df.tst$y, prob_noavg_mean[,1],lwd=1.5)
lines(df.tst$y, prob_noavg_med[,1], col="red",lwd=1.5)
lines(df.tst$y, prob_avg[,1], col="green",lwd=1.5)

legend(-7,0.14,c("Mean (no averaging)","Median (no averaging)","Bayes Average"),
       lty=c(1,1,1),col=c("black","red","green"),lwd=c(1.5,1.5,1.5), cex=0.7)


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
y_prob <- predict(rv, spambase.tst, type='prob')

# Check how well our predictions did by generating confusion matrix
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
