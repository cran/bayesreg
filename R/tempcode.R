# ============================================================================================================================
# Simple gradient descent fitting of generalized linear models
#[b, b0, L, gradb, gradb0] = br_FitGLMGD(X, y, model, xi, tau2, maxiter)
br.fit.GLM.gd <- function(df, model, xi, tau2, max.iter = 5e2)
{
  px = ncol(df) - 1
  py = nrow(df)
  theta = matrix(nrow = px+1, ncol = 1)
}