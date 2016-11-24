# ============================================================================================================================
# Declare the Bayesreg class
#setClass("bayesreg")
#library("BayesLogit")

#bayesreg <- function(formula, data, model='normal', prior='ridge', nsamples = 1e3, burnin = 1e3, thin = 5, display = F, groups = NA, tdof = 5, rankvars = T)
bayesreg <- function(formula, data, model='normal', prior='ridge', nsamples = 1e3, burnin = 1e3, thin = 5, tdof = 5)
{
  VERSION = '1.00'
  groups  = NA
  display = F
  rankvars = T
  
  printf <- function(...) cat(sprintf(...))
  
  # Return object
  rv = list()

  # -------------------------------------------------------------------    
  # modelibution types
  model = tolower(model)
  if (model == 'normal' || model == 'gaussian')
  {
    model = 'gaussian'
  }
  else if (model == 'laplace' || model == 'l1')
  {
    model = 'laplace'
  }
  else if (model == 'studentt' || model == 't')
  {
    model = 't'
  }
  else if (model == 'logistic' || model == "binomial")
  {
    model = 'logistic'
  }
  else
  {
    stop('Unknown target distribution \'',model,'\'')
  }
  
  # -------------------------------------------------------------------    
  # Process and set up the data from the model formula
  rv$terms <- terms(x = formula, data = data)
  
  mf = model.frame(formula = formula, data = data)
  rv$targetVar = names(mf)[1]
  if (model == 'logistic')
  {
    # Check to ensure target is a factor
    if (!is.factor(mf[,1]) || (!is.factor(mf[,1]) && length(levels(mf[,1])) != 2))
    {
      stop('Target variable must be a factor with two levels for logistic regression')
    }
    
    rv$ylevels <- levels(mf[,1])
    mf[,1] = as.numeric(mf[,1])
    mf[,1] = (as.numeric(mf[,1]) - 1)
  }
  
  # Otherwise check to ensure target is numeric
  else if (!is.numeric(mf[,1]))
  {
    stop('Target variable must not be a factor for linear regression; use model = "logistic" instead')
  }
  
  # If not logistic, check if targets have only 2 unique values -- if so, give a warning
  if (model != 'logistic' && length(unique(mf[,1])) == 2)
  {
    warning('Target variable takes on only two distinct values -- should this be a binary regression?')
  }
  
  y = mf[,1]
  X = as.matrix(model.matrix(formula, data=mf))
  X = X[,-1,drop=FALSE]
  
  # 
  n = nrow(X)
  p = ncol(X)

  # -------------------------------------------------------------------
  # Prior types
  if (prior == 'ridge' || prior == 'rr')
  {
    prior = 'rr'
  }
  else if (prior == 'horseshoe' || prior == 'hs')
  {
    prior = 'hs'
  }
  else if (prior == 'horseshoe+' || prior == 'hs+')
  {
    prior = 'hs+'
  }
  else if (prior != 'lasso')
  {
    stop('Unknown prior \'', prior, '\'')
  }
  
  # -------------------------------------------------------------------
  # Standardise data?
  stdX = bayesreg.standardise(X)
  X    = stdX$X
  XtX  = NA
  if (model == "gaussian" && p < 1e4)
  {
    XtX = t(X) %*% X
  }
  
  # Initial values
  ydiff   = 0
  b0      = 0
  b       = matrix(0,p,1)
  omega2  = matrix(1,n,1)
  sigma2  = 1

  tau2    = matrix(1,p,1)
  xi      = 1
  
  lambda2 = matrix(1,p,1)
  nu      = matrix(1,p,1)

  eta2    = matrix(1,p,1)
  phi     = matrix(1,p,1)

  kappa   = y - 1/2
  z       = y
  
  # Use Rue's MVN sampling algorithm
  mvnrue = T
  if (p/n >= 2)
  {
    # Switch to Bhatta's MVN sampling algorithm
    mvnrue = F
  }

  # Set up group shrinkage parameters if required
  ng = 0
  if (length(groups) > 1)
  {
    ngroups = length(unique(groups))
    delta2  = matrix(1,ngroups,1)
    chi     = matrix(1,ngroups,1)
  }
  else
  {
    ngroups = 1
    groups  = rep(1,1,p)
    delta2  = matrix(1,1,1)
    chi     = matrix(1,1,1)
  }

  # Setup return object
  rv$formula  = formula
  rv$model    = model
  rv$prior    = prior
  rv$nsamples = nsamples
  rv$burnin   = burnin
  rv$thin     = thin
  rv$groups   = groups
  rv$tdof     = tdof
  
  rv$n        = n
  rv$p        = p
  rv$tss      = sum((y - mean(y))^2)

  rv$beta0    = matrix(0,1,nsamples)
  rv$beta     = matrix(0,p,nsamples)
  rv$muBeta   = matrix(0,p,1)
  rv$muBeta0  = matrix(0,1,1)
  rv$sigma2   = matrix(0,1,nsamples)
  rv$muSigma2 = matrix(0,1,1)
  rv$tau2     = matrix(0,1,nsamples)
  #rv$xi       = matrix(0,1,nsamples)
  #rv$lambda2  = matrix(0,p,nsamples)
  #rv$nu       = matrix(0,p,nsamples)
  #rv$delta2   = matrix(0,ngroups,nsamples)
  #rv$chi      = matrix(0,ngroups,nsamples)
  
  rv$tstat    = matrix(0,p,1)
  
  class(rv)   = "bayesreg"
  
  # Print the banner  
  if (display)
  {
    printf('==========================================================================================\n');
    printf('|                   Bayesian Penalised Regression Estimation ver. %s                  |\n', VERSION);
    printf('|                         (c) Enes Makalic, Daniel F Schmidt. 2016                       |\n');
    printf('==========================================================================================\n');
  }

  # Main sampling loop
  k    = 0
  iter = 0
  while (k < nsamples)
  {
    # ===================================================================
    # Sample regression coefficients
    if (model == 'logistic')
    {
      z = kappa * omega2
    }
    
    # -------------------------------------------------------------------
    # Sample b0
    W       = sum(1 / omega2)
    muB0    = sum((z - X %*% b) / omega2) / W
    v       = sigma2 / W
    b0      = rnorm(1, muB0, sqrt(v))

    # -------------------------------------------------------------------
    # Sample regression coefficients
    bs  = bayesreg.sample_beta(X, z, mvnrue, b0, sigma2, tau2, lambda2 * delta2[groups], omega2, XtX)
    muB = bs$m
    b   = bs$x
    
    # ===================================================================
    # Sample the noise scale parameter sigma2
    muSigma2   = NA
    if (model != 'logistic')
    {
      e        = y - X %*% b - b0
      shape    = (n + p)/2
      scale    = sum( (e^2 / omega2)/2 ) + sum( (b^2 / delta2[groups] / lambda2 / tau2)/2 )
      sigma2   = 1/rgamma(1, shape=shape, scale=1/scale)
      
      muSigma2 = scale / (shape-1)
    }
    
    
    # ===================================================================
    # Sample 'omega's
    
    # -------------------------------------------------------------------
    # Laplace
    if (model == 'laplace')
    {
      e = y - X %*% b - b0
      mu = sqrt(2 * sigma2 / e^2)
      omega2 = 1 / bayesreg.rinvg(mu,1/2)
    }

    # -------------------------------------------------------------------
    # Student-t
    if (model == 't')
    {
      e = y - X %*% b - b0
      shape = (tdof+1)/2
      scale = (e^2/sigma2+tdof)/2
      omega2 = 1 / rgamma(n, shape=shape, scale=1/scale)
    }

    # -------------------------------------------------------------------
    # Logistic
    if (model == 'logistic')
    {
      omega2 = 1 / rpg.devroye(num = n, n = 1, z = b0 + X %*% b)
    }

    # ===================================================================
    # Sample the global shrinkage parameter tau2 (and L.V. xi)
    shape = (p+1)/2
    scale = 1/xi + sum(b^2 / lambda2 / delta2[groups]) / 2 / sigma2
    tau2  = 1 / rgamma(1, shape=shape, scale=1/scale)

    # Sample xi
    scale = 1 + 1/tau2
    xi    = 1 / rgamma(1, shape=1, scale=1/scale)
    
    
    # ===================================================================
    # Sample the lambda2's/nu's (if horseshoe, horseshoe+, horseshoe-grouped)
    if (prior == 'hs' || prior == 'hs+' || prior == 'hsge')
    {
      # Sample lambda2
      scale   = 1/nu + b^2 / 2 / tau2 / sigma2 / delta2[groups]
      lambda2 = 1 / rgamma(p, shape=1, scale=1/scale)
      
      # Sample nu -- horseshoe
      if (prior == 'hs' || prior == 'hsge')
      {
        scale = 1 + 1/lambda2
        nu    = 1 / rgamma(p, shape=1, scale=1/scale)
      }
      
      # Horseshoe+
      else if (prior == 'hs+')
      {
        # Sample nu
        scale = 1/eta2 + 1/lambda2
        nu    = 1 / rgamma(p, shape=1, scale=1/scale)
        
        # Sample eta2
        scale = 1/nu + 1/phi
        eta2  = 1 / rgamma(p, shape=1, scale=1/scale)
        
        # Sample phi
        scale = 1 + 1/eta2
        phi   = 1 / rgamma(p, shape=1, scale=1/scale)
      }
    }
    
    
    # ===================================================================
    # Sample the lambda2's (if lasso)
    if (prior == 'lasso')
    {
      mu      = sqrt(2 * sigma2 * tau2 / b^2)
      lambda2 = 1 / bayesreg.rinvg(mu, 1/2)
    }
    
    
    # ===================================================================
    # Sample the delta2's (if grouped horseshoe)
    if (prior == 'hsge' || prior == 'hsg')
    {
      # Sample delta2's
      for (i in 1:ngroups)
      {
        # Only sample delta2 for this group, if the group size is > 1
        ng = sum(groups == i)
        if (ng > 1)
        {
          # Sample delta2
          shape = (ng+1)/2
          scale = 1 / chi[i] + sum(b[groups==i]^2 / lambda2[groups==i]) / 2 / sigma2 / tau2
          delta2[i] = 1 / rgamma(1, shape=shape, scale=1/scale)
          
          # Sample chi
          scale  = 1 + 1/delta2[i]
          chi[i] = 1 / rgamma(1, shape=1, scale=1/scale)
        }
      }
    }
    
    # ===================================================================
    # Store the samples
    iter = iter+1
    if (iter > burnin)
    {
      # thinning
      if (!(iter%%thin))
      {
        k = k+1
        rv$beta0[k]    = b0
        rv$beta[,k]    = b
        rv$muBeta      = rv$muBeta + muB
        rv$muBeta0     = rv$muBeta0 + muB0
        rv$sigma2[k]   = sigma2
        rv$muSigma2    = rv$muSigma2 + muSigma2
        rv$tau2[k]     = tau2
        #rv$xi[k]       = xi
        #rv$lambda2[,k] = lambda2
        #rv$nu[,k]      = nu
        #rv$delta2[,k]  = delta2
        #rv$chi[,k]     = chi
      }
    }
  }
  
  #
  rv$muBeta   = rv$muBeta / nsamples
  rv$muBeta0  = rv$muBeta0 / nsamples
  rv$muSigma2 = rv$muSigma2 / nsamples
  
  # ===================================================================
  # Rank features, if requested
  if (rankvars)
  {
    # Run the BFR
    ranks = bayesreg.bfr(rv)

    # Determine the 75th percentile
    rv$varranks = rep(NA,p+1)
    q = apply(ranks,1,function(x) quantile(x,0.75))
    O = sort(q,index.return = T)
    O = O$ix
    
    j = 1
    k = 1
    for (i in 1:p)
    {
      if (i >= 2)
      {
        if (q[O[i]] != q[O[i-1]])
        {
          j = j+k
          k = 1
        }
        else
        {
          k = k+1
        }
      }
      rv$varranks[O[i]] = j
    }
  }
  else
  {
    rv$varranks = rep(NA, p+1)
  }
  
  # ===================================================================
  # Compute the t-statistics
  for (i in 1:p)
  {
    rv$tstat[i] = rv$muBeta[i] / sd(rv$beta[i,])
  }

  # ===================================================================
  # Compute DIC score and log-likelihoods
  rv$dic   = bayesreg.computedic(model, X, y, rv$beta, rv$beta0, tdof, rv);    
  if (model == 'logistic')
  {
    rv$logl  = -bayesreg.logregnlike(X, y, rv$muBeta, rv$muBeta0)
    rv$logl0 = -bayesreg.logregnlike(X, y, matrix(0,p,1), log(sum(y)/(n-sum(y))))
  }
  else
  {
    rv$logl  = -bayesreg.linregnlike(model, as.matrix(X), as.matrix(y), rv$muBeta, rv$muBeta0, rv$muSigma2, rv$tdof)
  }

  # ===================================================================
  # Rescale the coefficients
  if (p == 1)
  {
    rv$beta  <- t(as.matrix(apply(t(rv$beta), 1, function(x)(x / stdX$stdX))))
  }
  else
  {
    rv$beta  <- as.matrix(apply(t(rv$beta), 1, function(x)(x / stdX$stdX)))
  }
  rv$beta0 <- rv$beta0 - stdX$meanX %*% rv$beta
  
  rv$muBeta  <- rv$muBeta / t(stdX$stdX)
  rv$muBeta0 <- rv$muBeta0 - stdX$meanX %*% rv$muBeta

  rv$stdX = stdX
  
  # ===================================================================
  # Compute effective sample sizes
  rv$essfrac = rep(NA, p)
  for (i in 1:p)
  {
    e = bayesreg.ess(rv$beta[i,])
    rv$essfrac[i] = e$ESSfrac
  }
  rv$essfrac[ rv$essfrac > 100 ] = 100

  return(rv)
}


# ============================================================================================================================
# Sample the intercept
bayesreg.sample_beta0 <- function(X, z, b, sigma2, omega2)
{
  rv      = list()
  W       = sum(1 / omega2)
  rv$muB0 = sum((z - X %*% b) / omega2) / W
  v       = sigma2 / W
  
  rv$b0  = rnorm(1, rv$muB0, sqrt(v))
  
  return(rv)
}


# ============================================================================================================================
#
bayesreg.ess <- function(x)
{
  n = length(x)
  s = min(c(n-1, 2000))
  
  g = acf(x, s, plot=F)
  G = as.vector(g$acf[2:(s-1)]) + as.vector(g$acf[3:s])
  G = G < 0
  for (i in 1:length(G))
  {
    if (G[i])
    {
      break
    }
  }
  k = i
  if (k >= 2)
    V = g$acf[1] + 2 * sum(g$acf[2:i])
  else
    V = g$acf[1]
  
  ACT = V / g$acf[1]
  rv = list()
  rv$ESS = n/ ACT
  rv$ESSfrac = rv$ESS / n
  
  return(rv)
}


# ============================================================================================================================
# Bayesian Feature Ranking
bayesreg.bfr <- function(rv)
{
  ranks = matrix(0, rv$p, rv$nsamples)
  
  for (i in 1:rv$nsamples)
  {
    r = sort(-abs(rv$beta[,i]), index.return = T)
    ranks[r$ix,i] = 1:rv$p
  }
  
  return(ranks)
}   


# ============================================================================================================================
# predict function for bayesreg models
predict.bayesreg <- function(object, df, type = "linpred", bayesavg = FALSE, sum.stat = "mean", ...)
{
  # Error checking 
  if (sum.stat != "mean" && sum.stat != "median")
  {
    stop("The summary statistic must be either 'mean' or 'median'.")
  }
  if (type == "class" && object$model != "logistic")
  {
    stop("Class predictions are only available for logistic regressions.")
  }
  if (type != "linpred" && type != "prob" && type != "class")
  {
    stop("Type of prediction must be one of 'linpred', 'prob' or 'class'.")
  }
  
  # Build the fully specified formula using the covariates that were fitted
  f <- as.formula(paste("~",paste(attr(object$terms,"term.labels"),collapse="+")))
  
  # Extract the design matrix
  X = model.matrix(f, data=df)
  X = as.matrix(X[,-1])
  n = nrow(X)
  p = ncol(X)

  # Get y-data if it has been passed and is not NA
  if (!any(names(df) == object$targetVar) && type == "prob" && object$model != "logistic")
  {
    stop("You must provide a column of targets called '", object$targetVar, "' to predict probabilities for continuous distributions.")
  }
  
  if (any(names(df) == object$targetVar) && type == "prob")
  {
    y <- as.matrix(df[object$targetVar])
    if (any(is.na(y)) && type == "prob" && object$model != "logistic")
    {
      stop("Missing values in the target variable '", object$targetVar, "' not allowed when predicting probabilities for continuous distributions.")
    }
  }
  else {
    y <- NA
  }

  # Compute the linear predictor
  if (bayesavg == F)
  {
    #browser()
    if (sum.stat == "median")
    {
      medBeta  = apply(object$beta,1,function(x) quantile(x,c(0.5)))
      medBeta0 = quantile(object$beta0,0.5)
      yhat = X %*% as.vector(medBeta) + as.numeric(medBeta0)
    }
    else
    {
      yhat = X %*% as.vector(object$muBeta) + as.numeric(object$muBeta0)
    }
  }
  else
  {
    yhat = X %*% as.matrix(object$beta)
    yhat = t(as.matrix(apply(yhat,1,function(x) x + object$beta0)))
  }
  
  # Logistic reg -- class labels
  if (object$model == "logistic" && type == "class")
  {
    yhat = rowMeans(yhat)
    yhat[yhat>0] <- 2
    yhat[yhat<=0] <- 1
    yhat = factor(yhat, levels=c(1,2), labels=object$ylevels)

  # Logistic reg -- response prob
  } else if (object$model == "logistic" && type == "prob")
  {
    eps = 1e-16
    yhat = (1 / (1+exp(-yhat)) + eps)/(1+2*eps)
  
    # If binary and type="prob", and y has been passed, compute probabilities for the specific y values
    if (!any(is.na(y)) && type == "prob")
    {
      yhat[y == object$ylevels[1],] <- 1 - yhat[y == object$ylevels[1],]
    }
  }
  
  # Linear regression -- probability
  else if (object$model != "logistic" && type == "prob" && !any(is.na(y)))
  {
    # Gaussian probabilities
    if (object$model == "gaussian")
    {
      if (bayesavg == T)
      {
        yhat <- (as.matrix(apply(yhat, 2, function(x) (x - y)^2 )))
        yhat <- t(as.matrix(apply(yhat, 1, function(x) (1/2)*log(2*pi*object$sigma2) + x/2/object$sigma2)))
      }
      else
      {
        scale = as.numeric(object$muSigma2)
        yhat <- (1/2)*log(2*pi*scale) + (yhat - y)^2/2/scale
      }
    }

    # Laplace probabilities
    else if (object$model == "laplace")
    {
      if (bayesavg == T)
      {
        yhat <- as.matrix(apply(yhat, 2, function(x) abs(x - y) ))
        scale <- as.matrix(sqrt(object$sigma2/2))
        yhat <- t(as.matrix(apply(yhat, 1, function(x) log(2*scale) + x/scale)))
      }
      else
      {
        scale <- sqrt(as.numeric(object$muSigma2)/2)
        yhat <- log(2*scale) + abs(yhat - y)/scale
      }
    }
    
    # Student-t probabilities
    else
    {
      nu = object$tdof;
      if (bayesavg == T)
      {
        yhat <- (as.matrix(apply(yhat, 2, function(x) x - y)))^2
        scale = as.matrix(object$sigma2)
        yhat <- t(as.matrix(apply(yhat, 1, function(x) (-lgamma((nu+1)/2) + lgamma(nu/2) + log(pi*nu*scale)/2) + (nu+1)/2*log(1 + 1/nu*x/scale) )))
      }
      else
      {
        yhat <- (-lgamma((nu+1)/2) + lgamma(nu/2) + log(pi*nu*as.numeric(object$muSigma2))/2) + (nu+1)/2*log(1 + 1/nu*(yhat-y)^2/as.numeric(object$muSigma2))
      }
    }
    
    yhat <- exp(-yhat)
  }

  # If not class labels and averaging, also compute SE's and CI's
  if (type != "class" && bayesavg == T)
  {
    # Standard errors
    se = apply(yhat,1,sd)
    
    # 95% CIs
    CI = apply(yhat,1,function(x) quantile(x,c(0.05,0.5,0.95)))
    
    # Store results
    r = matrix(0, n, 4)
    if (sum.stat == "median" && type != "prob")
    {
      r[,1] = as.matrix(CI[2,])
    }
    else
    {
      r[,1] = as.matrix(rowMeans(yhat))
    }
    r[,2] = as.matrix(se)
    r[,3] = as.matrix(CI[1,])
    r[,4] = as.matrix(CI[3,])
    colnames(r) <- c("pred","se(pred)","CI 5%","CI 95%")
  }
  
  # Otherwise just return the class labels
  else
  {
    r = yhat
    if (object$model != "logistic")
    {
      colnames(r) <- "pred"
    }
  }
  
  return(r)
}


# ============================================================================================================================
# summary function for bayesreg models
summary.bayesreg <- function(object, sortrank = FALSE, displayor = FALSE, ...)
  {
    printf <- function(...) cat(sprintf(...))
    repchar <- function(c, n) paste(rep(c,n),collapse="")

    n   = object$n
    px  = object$p
    
    PerSD = F

    rv = list()

    # Error checking
    if (displayor == T && object$model != "logistic")
    {
      stop("Can only display cross-sectional odds-ratios for logistic models.")
    }
        
    # ===================================================================
    # Table symbols
    chline = '-'
    cvline = '|'
    cTT    = '+'
    cplus  = '+'
    bTT    = '+'
    
    # ===================================================================
    # Find length of longest variable name
    varnames = c(labels(object$muBeta)[[1]], "_cons")
    maxlen   = 12
    nvars    = length(varnames)
    for (i in 1:nvars)
    {
      if (nchar(varnames[i]) > maxlen)
      {
        maxlen = nchar(varnames[i])
      }
    }
    
    # ===================================================================
    # Display pre-table stuff
    #printf('%s%s%s\n', repchar('=', maxlen+1), '=', repchar('=', 76))
    printf('\n')
    
    if (object$model != 'logistic')
    {
      model = 'linear'
    } else {
      model = 'logistic'
    }
    
    if (object$model == "gaussian")
    {
      distr = "Gaussian"
    } else if (object$model == "laplace")
    {
      distr = "Laplace"
    }
    else if (object$model == "t")
    {
      distr = "Student-t"
    }
    else if (object$model == "logistic")
    {
      distr = "logistic"
    }
    
    if (object$prior == 'rr' || object$prior == 'ridge') {
      prior = 'ridge'
    } else if (object$prior == 'lasso') {      
      prior = 'lasso'
    } else if (object$prior == 'hs' || object$prior == 'horseshoe') {
      prior = 'horseshoe'
    } else if (object$prior == 'hs+' || object$prior == 'horseshoe+') {
      prior = 'horseshoe+'
    }

    str = sprintf('Bayesian %s %s regression', distr, prior)
    printf('%-64sNumber of obs   = %8d\n', str, n)
    printf('%-64sNumber of vars  = %8.0f\n', '', px);
    
    if (model == 'linear')
    {
      s2 = mean(object$sigma2)
      if (object$model == 't' && object$tdof > 2)
      {
        s2 = object$tdof/(object$tdof - 2) * s2
      }
      rsqr = 1 - (s2*(n-1)) / object$tss
      
      str = sprintf('MCMC Samples   = %6.0f', object$ns);
      printf('%-64sRoot MSE        = %8.5g\n', str, sqrt(s2));
      str = sprintf('MCMC Burnin    = %6.0f', object$burnin);
      printf('%-64sR-squared       = %8.4f\n', str, rsqr);    
      str = sprintf('MCMC Thinning  = %6.0f', object$thin);
      printf('%-64sDIC             = %8.5g\n', str, object$dic);
      
      rv$rootmse = sqrt(s2)
      rv$rsqr    = rsqr
      rv$dic     = object$dic
      rv$logl    = object$logl
    }
    else if (model == 'logistic')
    {
      logl = object$logl
      logl0 = object$logl0
      r2 = 1 - logl / logl0
      
      str = sprintf('MCMC Samples   = %6.0f', object$ns);
      printf('%-64sLog. Likelihood = %8.5g\n', str, object$logl);
      str = sprintf('MCMC Burnin    = %6.0f', object$burnin);
      printf('%-64sPseudo R2       = %8.4f\n', str, r2);    
      str = sprintf('MCMC Thinning  = %6.0f', object$thin);
      printf('%-64sDIC             = %8.5g\n', str, object$dic);
      
      rv$logl   = logl
      rv$p_rsqr = r2
      rv$dic    = object$dic
    }
    printf('\n')
    
    # ===================================================================
    # Table Header
    fmtstr = sprintf('%%%ds', maxlen);
    printf('%s%s%s\n', repchar(chline, maxlen+1), cTT, repchar(chline, 76))
    tmpstr = sprintf(fmtstr, 'Parameter');
    if (model == 'linear' || displayor == F)
    {
      printf('%s %s  %10s %10s    [95%% Cred. Interval] %10s %7s %9s\n', tmpstr, cvline, 'mean(Coef)', 'std(Coef)', 'tStat', 'Rank', 'ESS');
    } else if (model == 'logistic' && displayor == T) {
      printf('%s %s  %10s %10s    [95%% Cred. Interval] %10s %7s %9s\n', tmpstr, cvline, 'median(OR)', 'std(OR)', 'tStat', 'Rank', 'ESS');
    }
    printf('%s%s%s\n', repchar(chline, maxlen+1), cTT, repchar(chline, 76));

    # ===================================================================
    if (PerSD) {
      beta0 = object$beta0
      beta  = object$beta
      
      object$beta0  = beta0 + object$stdX$meanX %*% object$beta
      object$beta   = apply(t(object$beta), 1, function(x)(x * object$stdX$stdX/sqrt(n)))    
      object$muBeta = object$muBeta * t(as.matrix(object$stdX$stdX)) / sqrt(n)
    }

    # ===================================================================
    # Return values
    rv$muCoef   = matrix(0,px+1,1, dimnames = list(varnames))
    rv$seCoef   = matrix(0,px+1,1, dimnames = list(varnames))
    rv$CICoef   = matrix(0,px+1,2, dimnames = list(varnames))
    
    if (model == "logistic")
    {
      rv$medOR  = matrix(0,px+1,1, dimnames = list(varnames))
      rv$seOR   = matrix(0,px+1,1, dimnames = list(varnames))
      rv$CIOR   = matrix(0,px+1,2, dimnames = list(varnames))
    }

    rv$tStat    = matrix(0,px+1,1, dimnames = list(varnames))
    rv$nStars   = matrix(0,px+1,1, dimnames = list(varnames))
    rv$ESS      = matrix(0,px+1,1, dimnames = list(varnames))
    rv$rank     = matrix(0,px+1,1, dimnames = list(varnames))

    # Variable information
    rv$rank[1:(px+1)] = as.matrix(object$varranks)

    # If sorting by BFR ranks
    if (sortrank == T)
    {
      O = sort(object$varranks, index.return = T)
      indices = O$ix
      indices[px+1] = px+1
    }    
    # Else sorted by the order they were passed in
    else {
      indices = (1:(px+1))
    }
    
    for (i in 1:(px+1))
    {
      k = indices[i]
      kappa = NA
      
      # Regression variable
      if (k <= px)
      {
        s = object$beta[k,]
        mu = object$muBeta[k]

        # Calculate shrinkage proportion, if possible
        kappa = object$tstat[k]
      }
      
      # Intercept
      else if (k == (px+1))
      {
        s = object$beta0
        mu = mean(s)
      }
      
      # Compute credible intervals/standard errors for beta's
      std_err = sd(s)
      qlin = quantile(s, c(0.025, 0.25, 0.75, 0.975))
      qlog = quantile(exp(s), c(0.025, 0.25, 0.75, 0.975))
      q = qlin

      rv$muCoef[k] = mu
      rv$seCoef[k] = std_err
      rv$CICoef[k,] = c(qlin[1],qlin[4])

      # Compute credible intervals/standard errors for OR's
      if (model == 'logistic')
      {
        med_OR = median(exp(s))
        std_err_OR = (qlog[4]-qlog[1])/2/1.96
        
        rv$medOR[k] = med_OR
        rv$seOR[k]  = std_err_OR
        rv$CIOR[k,] = c(qlog[1],qlog[4])

        # If display ORs, use these instead
        if (displayor)
        {
          mu = med_OR
          std_err = std_err_OR
          q = qlog
        }
      }
      rv$tStat[k] = kappa
      
      # Display results
      tmpstr = sprintf(fmtstr, varnames[k])
      if (is.na(kappa))
        tstat = '         .'
      else tstat = sprintf('%10.3f', kappa)
      
      if (is.na(object$varranks[k]))
        rank = '      .'
      else
        rank = sprintf('%7d', object$varranks[k])

      printf('%s %s %11.5f %10.5f   %10.5f %10.5f %s %s ', tmpstr, cvline, mu, std_err, q[1], q[4], tstat, rank);

      # Model selection scores
      # Test if 75% CI includes 0
      if ( k <= px && ( (qlin[2] > 0 && qlin[3] > 0) || (qlin[2] < 0 && qlin[3] < 0) ) )    
      {
        printf('*')
        rv$nStars[k] = rv$nStars[k] + 1
      }
      else 
        printf(' ')
      
      # Test if 95% CI includes 0
      if ( k <= px && ( (qlin[1] > 0 && qlin[4] > 0) || (qlin[1] < 0 && qlin[4] < 0) ) )    
      {
        printf('*')
        rv$nStars[k] = rv$nStars[k] + 1
      }
      else
        printf(' ')

      # Display ESS-frac
      if(k > px)
        printf('%7s', '.')
      else
        printf('%7.1f', object$essfrac[k]*100)

      rv$ESS[k] = object$essfrac[k]*100
      
      printf('\n');
    }
    
    printf('%s%s%s\n\n', repchar(chline, maxlen+1), cTT, repchar(chline, 76));
    
    return(rv)
}


# ============================================================================================================================
# function to standardise columns of X to have mean zero and unit length
bayesreg.standardise <- function(X)
{
  n = nrow(X)
  p = ncol(X)
  
  # 
  r       = list()
  r$X     = X
  if (p > 1)
  {
    r$meanX = colMeans(X)
  } else
  {
    r$meanX = mean(X)
  }
  r$stdX  = t(apply(X,2,sd)) * sqrt(n-1)
  
  # Perform the standardisation
  if (p == 1)
  {
    r$X <- as.matrix(apply(X,1,function(x)(x - r$meanX)))
    r$X <- as.matrix(apply(r$X,1,function(x)(x / r$stdX)))
  } else
  {
    r$X <- t(as.matrix(apply(X,1,function(x)(x - r$meanX))))
    r$X <- t(as.matrix(apply(r$X,1,function(x)(x / r$stdX))))
  }
  
  return(r)
}


# ============================================================================================================================
# Sample the regression coefficients
bayesreg.sample_beta <- function(X, z, mvnrue, b0, sigma2, tau2, lambda2, omega2, XtX)
{
  alpha  = (z - b0)
  Lambda = sigma2 * tau2 * lambda2
  sigma  = sqrt(sigma2)
  
  # Use Rue's algorithm
  if (mvnrue)
  {
    # If XtX is not precomputed
    if (any(is.na(XtX)))
    {
      omega = sqrt(omega2)
      X0    = apply(X,2,function(x)(x/omega))
      bs    = bayesreg.fastmvg2_rue(X0/sigma, alpha/sigma/omega, Lambda)
    }
    
    # XtX is precomputed (Gaussian only)
    else {
      bs    = bayesreg.fastmvg2_rue(X/sigma, alpha/sigma, Lambda, XtX/sigma2)
    }
  }
  
  # Else use Bhat. algorithm
  else
  {
    omega = sqrt(omega2)
    X0    = apply(X,2,function(x)(x/omega))
    bs    = bayesreg.fastmvg_bhat(X0/sigma, alpha/sigma/omega, Lambda)
  }
  
  return(bs)
}


# ============================================================================================================================
# function to generate multivariate normal random variates using Rue's algorithm
bayesreg.fastmvg2_rue <- function(Phi, alpha, d, PtP = NA)
{
  Phi   = as.matrix(Phi)
  alpha = as.matrix(alpha)
  r     = list()
  
  # If PtP not precomputed
  if (any(is.na(PtP)))
  {
    PtP = t(Phi) %*% Phi
  }
  
  p     = ncol(Phi)
  if (length(d) > 1)
  {
    Dinv  = diag(as.vector(1/d))
  }
  else
  {
    Dinv   = 1/d
  }
  L     = t(chol(PtP + Dinv))
  v     = forwardsolve(L, t(Phi) %*% alpha)
  r$m   = backsolve(t(L), v)
  w     = backsolve(t(L), rnorm(p,0,1))
  
  r$x   = r$m + w
  return(r)
}


# ============================================================================================================================
# function to generate multivariate normal random variates using Bhat. algorithm
bayesreg.fastmvg_bhat <- function(Phi, alpha, d)
{
  d     = as.matrix(d)
  p     = ncol(Phi)
  n     = nrow(Phi)
  r     = list()
  
  u     = as.matrix(rnorm(p,0,1)) * sqrt(d)
  delta = as.matrix(rnorm(n,0,1))
  v     = Phi %*% u + delta
  Dpt   = (apply(Phi, 1, function(x)(x*d)))
  W     = Phi %*% Dpt + diag(1,n)
  w     = solve(W,(alpha-v))
  
  r$x   = u + Dpt %*% w
  r$m   = r$x

  return(r)
}


# ============================================================================================================================
# rinvg
bayesreg.rinvg <- function(mu, lambda)
{
  lambda = 1/lambda
  p = length(mu)
  
  V = rnorm(p,mean=0,sd=1)^2
  out = mu + 1/2*mu/lambda * (mu*V - sqrt(4*mu*lambda*V + mu^2*V^2))
  out[out<1e-16] = 1e-16
  z = runif(p)
  
  l = z >= mu/(mu+out)
  out[l] = mu[l]^2 / out[l]
  
  return(out)
}


# ============================================================================================================================
# compute the DIC score for the model
bayesreg.computedic <- function(error, X, y, beta, beta0, tdof, rv)
{
  nsamples = ncol(beta)
  
  dic = 0
  negll = 0
  for (i in 1:nsamples)
  {
    if (error != 'logistic')
    {
      negll = bayesreg.linregnlike(error, as.matrix(X), as.matrix(y), beta[,i], beta0[i], rv$sigma2[i], tdof)
    }
    else
    {
      negll = bayesreg.logregnlike(as.matrix(X), as.matrix(y), beta[,i], beta0[i])
    }
    dic = dic + negll
  }
  dic = dic / nsamples
  
  if (error != 'logistic')
  {
    dic = -2*dic + bayesreg.linregnlike(error, as.matrix(X), as.matrix(y), rv$muBeta, rv$muBeta0, rv$muSigma2, tdof)
  }
  else
  {
    dic = -2*dic + bayesreg.logregnlike(as.matrix(X), as.matrix(y), rv$muBeta, rv$muBeta0)
  }
  
  return(dic)
}

# ============================================================================================================================
# compute the negative log-likelihood
bayesreg.linregnlike <- function(error, X, y, b, b0, s2, tdof = NA)
{
  n = nrow(y)

  y = as.matrix(y)
  X = as.matrix(X)
  b = as.matrix(b)
  b0 = as.numeric(b0)
  
  e = (y - X %*% as.matrix(b) - as.numeric(b0))
  
  if (error == 'gaussian')
  {
    negll = n/2*log(2*pi*s2) + 1/2/s2*t(e) %*% e
  }
  else if (error == 'laplace')
  {
    scale <- sqrt(as.numeric(s2)/2)
    negll = n*log(2*scale) + sum(abs(e))/scale
  }
  else if (error == 't')
  {
    nu = tdof;
    negll = n*(-lgamma((nu+1)/2) + lgamma(nu/2) + log(pi*nu*s2)/2) + (nu+1)/2*sum(log(1 + 1/nu*e^2/as.numeric(s2)))
  }
  
  return(negll)
}

# ============================================================================================================================
# compute the negative log-likelihood for logistic regression
bayesreg.logregnlike <- function(X, y, b, b0)
{
  y = as.matrix(y)
  X = as.matrix(X)
  b = as.matrix(b)
  b0 = as.numeric(b0)
  
  eps = exp(-36)
  lowerBnd = log(eps)
  upperBnd = -lowerBnd
  muLims = c(eps, 1-eps)
  
  eta = as.numeric(b0) + X %*% as.matrix(b)
  eta[eta < lowerBnd] = lowerBnd
  eta[eta > upperBnd] = upperBnd
  
  mu = 1 / (1 + exp(-eta))
  
  mu[mu < eps] = eps
  mu[mu > (1-eps)] = (1-eps)
  
  negll = -sum(y*log(mu)) -
    sum((1.0 - y)*log(1.0 - mu))
  
  return(negll)
}