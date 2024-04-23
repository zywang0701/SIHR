#' Inference for linear combination of the regression vector in high dimensional
#' generalized linear regression
#'
#' @param X Design matrix, of dimension \eqn{n} x \eqn{p}
#' @param y Outcome vector, of length \eqn{n}
#' @param loading.mat Loading matrix, nrow=\eqn{p}, each column corresponds to a
#'   loading of interest
#' @param model The high dimensional regression model, either \code{"linear"} or
#'   \code{"logistic"} or \code{"logistic_alter"}
#' @param intercept Should intercept be fitted for the initial estimator
#'   (default = \code{TRUE})
#' @param intercept.loading Should intercept term be included for the loading
#'   (default = \code{FALSE})
#' @param beta.init The initial estimator of the regression vector (default =
#'   \code{NULL})
#' @param lambda The tuning parameter in fitting initial model. If \code{NULL},
#'   it will be picked by cross-validation. (default = \code{NULL})
#' @param mu The dual tuning parameter used in the construction of the
#'   projection direction. If \code{NULL} it will be searched automatically.
#'   (default = \code{NULL})
#' @param prob.filter The threshold of estimated probabilities for filtering
#'   observations in logistic regression. (default = 0.05)
#' @param rescale The factor to enlarge the standard error to account for the
#'   finite sample bias. (default = 1.1)
#' @param verbose Should intermediate message(s) be printed. (default =
#'   \code{FALSE})
#'
#' @return
#' \item{est.plugin.vec}{The vector of plugin(biased) estimators for the
#' linear combination of regression coefficients, length of
#' \code{ncol(loading.mat)}; each corresponding to a loading of interest}
#' \item{est.debias.vec}{The vector of bias-corrected estimators for the linear
#' combination of regression coefficients, length of \code{ncol(loading.mat)};
#' each corresponding to a loading of interest}
#' \item{se.vec}{The vector of standard errors of the bias-corrected estimators,
#' length of \code{ncol(loading.mat)}; each corresponding to a loading of interest}
#' \item{proj.mat}{The matrix of projection directions; each column corresponding
#' to a loading of interest.}
#'
#' @export
#' @import CVXR glmnet
#' @importFrom stats coef dnorm median pnorm qnorm symnum
#'
#' @examples
#' X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
#' y <- -0.5 + X[, 1] * 0.5 + X[, 2] * 1 + rnorm(100)
#' loading1 <- c(1, 1, rep(0, 3))
#' loading2 <- c(-0.5, -1, rep(0, 3))
#' loading.mat <- cbind(loading1, loading2)
#' Est <- LF(X, y, loading.mat, model = "linear")
#'
#' ## compute confidence intervals
#' ci(Est, alpha = 0.05, alternative = "two.sided")
#'
#' ## summary statistics
#' summary(Est)
LF <- function(X, y, loading.mat, model = c("linear", "logistic", "logistic_alter"),
               intercept = TRUE, intercept.loading = FALSE, beta.init = NULL, lambda = NULL,
               mu = NULL, prob.filter = 0.05, rescale = 1.1, verbose = FALSE) {
  ## Check arguments
  model <- match.arg(model)
  X <- as.matrix(X)
  y <- as.vector(y)
  loading.mat <- as.matrix(loading.mat)
  check.args.LF(
    X = X, y = y, loading.mat = loading.mat, model = model, intercept = intercept,
    intercept.loading = intercept.loading, beta.init = beta.init, lambda = lambda,
    mu = mu, rescale = rescale, prob.filter = prob.filter,
    verbose = verbose
  )
  loading_include_intercept <- (ncol(X) == (nrow(loading.mat) - 1))


  # Preparation -------------------------------------------------------------
  ## Centralize X
  X_means <- colMeans(X)
  X <- scale(X, center = TRUE, scale = F)

  ## Specify relevant functions
  funs.all <- relevant.funs(intercept = intercept, model = model)
  train.fun <- funs.all$train.fun
  pred.fun <- funs.all$pred.fun
  deriv.fun <- funs.all$deriv.fun
  weight.fun <- funs.all$weight.fun
  cond_var.fun <- funs.all$cond_var.fun

  ## Initial lasso estimator of beta
  if (is.null(beta.init)) beta.init <- train.fun(X, y, lambda = lambda)$lasso.est
  beta.init <- as.vector(beta.init)
  sparsity <- sum(abs(beta.init) > 1e-4)

  ## Obtain relevant values
  if (intercept) X <- cbind(1, X)
  n <- nrow(X)
  p <- ncol(X)
  pred <- as.vector(pred.fun(X %*% beta.init))
  deriv <- as.vector(deriv.fun(X %*% beta.init))
  weight <- as.vector(weight.fun(X %*% beta.init))
  cond_var <- as.vector(cond_var.fun(pred, y, sparsity))

  # Filtering Observations --------------------------------------------------
  ## For logistic regressions, we should filter out those observations whose
  ## values of estimated probability are too extreme.
  idx <- rep(TRUE, n)
  if (model != "linear") {
    idx <- as.logical((pred > prob.filter) * (pred < (1 - prob.filter)))
    if (mean(idx) < 0.8) message("More than 20 % observations are filtered out
                                 as their estimated probabilities are too close to the
                                 boundary 0 or 1.")
  }
  X.filter <- X[idx, , drop = F]
  y.filter <- y[idx]
  weight.filter <- weight[idx]
  deriv.filter <- deriv[idx]
  n.filter <- nrow(X.filter)

  # Establish Inference for loadings ----------------------------------------
  n.loading <- ncol(loading.mat)
  est.plugin.vec <- est.debias.vec <- se.vec <- rep(NA, n.loading)
  proj.mat <- matrix(NA, nrow = p, ncol = n.loading)

  for (i.loading in 1:n.loading) {
    if (verbose) cat(sprintf("---> Computing for loading (%i/%i)... \n", i.loading, n.loading))
    ## Adjust loading
    loading <- as.vector(loading.mat[, i.loading])
    if (loading_include_intercept) {
      loading <- loading
    } else {
      if (intercept) {
        if (intercept.loading) {
          loading <- loading - X_means
          loading <- c(1, loading)
        } else {
          loading <- c(0, loading)
        }
      }
    }
    loading.norm <- sqrt(sum(loading^2))

    ## Correction Direction
    direction <- compute_direction(loading, X.filter, weight.filter, deriv.filter, mu, verbose)

    ## Bias Correction
    est.plugin <- sum(beta.init * loading)
    correction <- mean((weight * (y - pred) * X) %*% direction)
    est.debias <- est.plugin + correction * loading.norm

    ## Compute SE and Construct CI
    V <- sum(((sqrt(weight^2 * cond_var) * X) %*% direction)^2) / n * loading.norm^2
    if ((n > 0.9 * p) & (model == "linear")) se <- sqrt(V / n) else se <- rescale * sqrt(V / n)

    ## Store Infos
    est.plugin.vec[i.loading] <- est.plugin
    est.debias.vec[i.loading] <- est.debias
    se.vec[i.loading] <- se
    proj.mat[, i.loading] <- direction * loading.norm
  }

  obj <- list(
    est.plugin.vec = est.plugin.vec,
    est.debias.vec = est.debias.vec,
    se.vec = se.vec,
    proj.mat = proj.mat
  )

  class(obj) <- "LF"
  obj
}
