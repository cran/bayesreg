#' This is a comprehensive, user-friendly package implementing the state-of-the-art in 
#' Bayesian linear regression and Bayesian logistic regression. Features of the toolbox include:
#' \itemize{
#'  \item Supports Gaussian, Laplace, Student-t and logistic binary data models.
#'  \item Efficient and numerically stable implementations of Bayesian ridge, Bayesian lasso, horseshoe and horseshoe+ regression.
#'  \item Provides variable ranking and importance, credible intervals and diagnostics such as the widely applicable information criterion.
#'  \item Factor variables are automatically grouped together and additional shrinkage is applied to the set of indicator variables to which they expand.
#'  \item Prediction tools for generating credible intervals and Bayesian averaging of predictions.
#' }
#' The lasso, horseshoe and horseshoe+ priors are recommended for data sets where the number of 
#' predictors is greater than the sample size. The non-Gaussian models are based on scale-mixture representations.
#' Logistic regression utilises the Polya-gamma sampler implemented in the \code{pgdraw} package.
#'
#' To support analysis of data with outliers, we provide two heavy-tailed error models in our 
#' implementation of Bayesian linear regression: Laplace and Student-t distribution errors. 
#' The widely applicable information criterion (WAIC) is routinely calculated and displayed
#' to assist users in selecting an appropriate prior distribution for their particular problem, i.e., choice of regularisation or data model.
#' Most features are straightforward to use.
#' 
#' Further information on the particular algorithms/methods implemented in this package provided
#' by the literature referenced below.
#' 
#' @title Getting started with the bayesreg package
#' @docType package
#' @author Daniel Schmidt \email{daniel.schmidt@@monash.edu} 
#' 
#' Faculty of Information Technology, Monash University, Australia
#'
#' Enes Makalic \email{emakalic@@unimelb.edu.au}
#' 
#' Centre for Epidemiology and Biostatistics, The University of Melbourne, Australia
#' 
#' @note     To cite this package please reference: 
#'
#' Makalic, E. & Schmidt, D. F.
#' High-Dimensional Bayesian Regularised Regression with the BayesReg Package
#' arXiv:1611.06649 [stat.CO], 2016 \url{https://arxiv.org/pdf/1611.06649.pdf}
#' 
#' A MATLAB-compatible implementation of this package can be obtained from:
#' 
#' \url{https://au.mathworks.com/matlabcentral/fileexchange/60823}
#'
#' @seealso \code{\link{bayesreg}}
#' @name bayesreg-package
#' @references 
#'
#'    Bhadra, A.; Datta, J.; Polson, N. G. & Willard, B. 
#'    The Horseshoe+ Estimator of Ultra-Sparse Signals 
#'    Bayesian Analysis, 2016
#'
#'    Bhattacharya, A.; Chakraborty, A. & Mallick, B. K. 
#'    Fast sampling with Gaussian scale-mixture priors in high-dimensional regression 
#'    arXiv:1506.04778, 2016
#'
#'    Carvalho, C. M.; Polson, N. G. & Scott, J. G. 
#'    The horseshoe estimator for sparse signals 
#'    Biometrika, Vol. 97, pp. 465-480, 2010
#'
#'    Makalic, E. & Schmidt, D. F. 
#'    A Simple Sampler for the Horseshoe Estimator 
#'    IEEE Signal Processing Letters, Vol. 23, pp. 179-182, 2016
#'
#'    Park, T. & Casella, G. 
#'    The Bayesian Lasso 
#'    Journal of the American Statistical Association, Vol. 103, pp. 681-686, 2008
#'
#'    Polson, N. G.; Scott, J. G. & Windle, J. 
#'    Bayesian inference for logistic models using Polya-Gamma latent variables 
#'    Journal of the American Statistical Association, Vol. 108, pp. 1339-1349, 2013
#'
#'    Rue, H. 
#'    Fast sampling of Gaussian Markov random fields 
#'    Journal of the Royal Statistical Society (Series B), Vol. 63, pp. 325-338, 2001
#'
#'    Xu, Z., Schmidt, D.F., Makalic, E., Qian, G. & Hopper, J.L.
#'    Bayesian Grouped Horseshoe Regression with Application to Additive Models
#'    AI 2016: Advances in Artificial Intelligence, pp. 229-240, 2016
#' 
#' @examples 
#' # -----------------------------------------------------------------
#' # Example 1: Gaussian regression
#' X = matrix(rnorm(100*20),100,20)
#' b = matrix(0,20,1)
#' b[1:5] = c(5,4,3,2,1)
#' y = X %*% b + rnorm(100, 0, 1)
#' 
#' df <- data.frame(X,y)
#' rv.lm <- lm(y~.,df)                        # Regular least-squares
#' summary(rv.lm)
#' 
#' rv.hs <- bayesreg(y~.,df,prior="hs")       # Horseshoe regression
#' rv.hs.s <- summary(rv.hs)
#' 
#' # Expected squared prediction error for least-squares
#' coef_ls = coef(rv.lm)
#' as.numeric(sum( (as.matrix(coef_ls[-1]) - b)^2 ) + coef_ls[1]^2)
#' 
#' # Expected squared prediction error for horseshoe
#' as.numeric(sum( (rv.hs$mu.beta - b)^2 ) + rv.hs$mu.beta0^2)
#' 
#' 
#' # -----------------------------------------------------------------
#' # Example 2: Gaussian v Student-t robust regression
#' X = 1:10;
#' y = c(-0.6867, 1.7258, 1.9117, 6.1832, 5.3636, 7.1139, 9.5668, 10.0593, 11.4044, 6.1677);
#' df = data.frame(X,y)
#' 
#' # Gaussian ridge
#' rv.G <- bayesreg(y~., df, model = "gaussian", prior = "ridge", n.samples = 1e3)
#' 
#' # Student-t ridge
#' rv.t <- bayesreg(y~., df, model = "t", prior = "ridge", t.dof = 5, n.samples = 1e3)
#' 
#' # Plot the different estimates with credible intervals
#' plot(df$X, df$y, xlab="x", ylab="y")
#' 
#' yhat_G <- predict(rv.G, df, bayes.avg=TRUE)
#' lines(df$X, yhat_G[,1], col="blue", lwd=2.5)
#' lines(df$X, yhat_G[,3], col="blue", lwd=1, lty="dashed")
#' lines(df$X, yhat_G[,4], col="blue", lwd=1, lty="dashed")
#' 
#' yhat_t <- predict(rv.t, df, bayes.avg=TRUE)
#' lines(df$X, yhat_t[,1], col="darkred", lwd=2.5)
#' lines(df$X, yhat_t[,3], col="darkred", lwd=1, lty="dashed")
#' lines(df$X, yhat_t[,4], col="darkred", lwd=1, lty="dashed")
#' 
#' legend(1,11,c("Gaussian","Student-t (dof=5)"),lty=c(1,1),col=c("blue","darkred"),
#'        lwd=c(2.5,2.5), cex=0.7)
#' 
#' \dontrun{
#' # -----------------------------------------------------------------
#' # Example 3: Logistic regression on spambase
#' data(spambase)
#'   
#' # bayesreg expects binary targets to be factors
#' spambase$is.spam <- factor(spambase$is.spam)
#' 
#' # First take a subset of the data (1/10th) for training, reserve the rest for testing
#' spambase.tr  = spambase[seq(1,nrow(spambase),10),]
#' spambase.tst = spambase[-seq(1,nrow(spambase),10),]
#'   
#' # Fit a model using logistic horseshoe for 2,000 samples
#' rv <- bayesreg(is.spam ~ ., spambase.tr, model = "logistic", prior = "horseshoe", n.samples = 2e3)
#'   
#' # Summarise, sorting variables by their ranking importance
#' rv.s <- summary(rv,sort.rank=TRUE)
#'   
#' # Make predictions about testing data -- get class predictions and class probabilities
#' y_pred <- predict(rv, spambase.tst, type='class')
#'   
#' # Check how well did our predictions did by generating confusion matrix
#' table(y_pred, spambase.tst$is.spam)
#'   
#' # Calculate logarithmic loss on test data
#' y_prob <- predict(rv, spambase.tst, type='prob')
#' cat('Neg Log-Like for no Bayes average, posterior mean estimates: ', sum(-log(y_prob[,1])), '\n')
#' y_prob <- predict(rv, spambase.tst, type='prob', sum.stat="median")
#' cat('Neg Log-Like for no Bayes average, posterior median estimates: ', sum(-log(y_prob[,1])), '\n')
#' y_prob <- predict(rv, spambase.tst, type='prob', bayes.avg=TRUE)
#' cat('Neg Log-Like for Bayes average: ', sum(-log(y_prob[,1])), '\n')
#' }
NULL