#' Inference for quadratic forms of the regression vector in high dimensional
#' generalized linear regressions
#'
#' @param X Design matrix, of dimension \eqn{n} x \eqn{p}
#' @param y Outcome vector, of length \eqn{n}
#' @param G The set of indices, \code{G} in the quadratic form
#' @param A The matrix A in the quadratic form, of dimension
#'   \eqn{|G|\times}\eqn{|G|}. If \code{NULL} A would be set as the
#'   \eqn{|G|\times}\eqn{|G|} submatrix of the population covariance matrix
#'   corresponding to the index set \code{G} (default = \code{NULL})
#' @param model The high dimensional regression model, either \code{"linear"} or
#'   \code{"logistic"} or \code{"logistic_alter"}
#' @param intercept Should intercept be fitted for the initial estimator
#'   (default = \code{TRUE})
#' @param beta.init The initial estimator of the regression vector (default =
#'   \code{NULL})
#' @param split Sampling splitting or not for computing the initial estimator.
#'   It take effects only when \code{beta.init =  NULL}. (default = \code{TRUE})
#' @param lambda The tuning parameter in fitting initial model. If \code{NULL},
#'   it will be picked by cross-validation. (default = \code{NULL})
#' @param mu The dual tuning parameter used in the construction of the
#'   projection direction. If \code{NULL} it will be searched automatically.
#'   (default = \code{NULL})
#' @param prob.filter The threshold of estimated probabilities for filtering
#'   observations in logistic regression. (default = 0.05)
#' @param rescale The factor to enlarge the standard error to account for the
#'   finite sample bias. (default = 1.1)
#' @param tau The enlargement factor for asymptotic variance of the
#'   bias-corrected estimator to handle super-efficiency. It allows for a scalar
#'   or vector. (default = \code{c(0.25,0.5,1)})
#' @param verbose Should intermediate message(s) be printed. (default =
#'   \code{FALSE})
#'
#' @return
#' \item{est.plugin}{The plugin(biased) estimator for the quadratic form of the
#' regression vector restricted to \code{G}}
#' \item{est.debias}{The bias-corrected estimator of the quadratic form of the
#' regression vector}
#' \item{se}{Standard errors of the bias-corrected estimator,
#' length of \code{tau}; corrsponding to different values of \code{tau}}
#' @export
#'
#' @import CVXR glmnet
#' @importFrom stats coef dnorm median pnorm qnorm symnum
#' @examples
#' X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
#' y <- X[, 1] * 0.5 + X[, 2] * 1 + rnorm(100)
#' G <- c(1, 2)
#' A <- matrix(c(1.5, 0.8, 0.8, 1.5), nrow = 2, ncol = 2)
#' Est <- QF(X, y, G, A, model = "linear")
#' ## compute confidence intervals
#' ci(Est, alpha = 0.05, alternative = "two.sided")
#'
#' ## summary statistics
#' summary(Est)
QF <- function(X, y, G, A = NULL, model = c("linear", "logistic", "logistic_alter"),
               intercept = TRUE, beta.init = NULL, split = TRUE, lambda = NULL, mu = NULL,
               prob.filter = 0.05, rescale = 1.1, tau = c(0.25, 0.5, 1), verbose = FALSE) {
  ## Check arguments
  model <- match.arg(model)
  X <- as.matrix(X)
  y <- as.vector(y)
  G <- as.vector(G)
  nullA <- ifelse(is.null(A), TRUE, FALSE)
  check.args.QF(
    X = X, y = y, G = G, A = A, model = model, intercept = intercept, beta.init = beta.init,
    split = split, lambda = lambda, mu = mu, prob.filter = prob.filter,
    rescale = rescale, tau = tau, verbose = verbose
  )

  # Preparation -------------------------------------------------------------
  ### Centralize X ###
  X <- scale(X, center = TRUE, scale = F)

  ## Specify relevant functions
  funs.all <- relevant.funs(intercept = intercept, model = model)
  train.fun <- funs.all$train.fun
  pred.fun <- funs.all$pred.fun
  deriv.fun <- funs.all$deriv.fun
  weight.fun <- funs.all$weight.fun
  cond_var.fun <- funs.all$cond_var.fun

  ## Initial lasso estimator of beta
  if (is.null(beta.init)) {
    if (split) {
      idx1 <- sample(1:nrow(X), size = round(nrow(X) / 2), replace = F)
      idx2 <- setdiff(1:nrow(X), idx1)
      X1 <- X[idx1, , drop = F]
      y1 <- y[idx1]
      X <- X[idx2, , drop = F]
      y <- y[idx2]
    } else {
      X1 <- X
      y1 <- y
    }
    beta.init <- train.fun(X1, y1, lambda = lambda)$lasso.est
  }
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
  ## Adjust A and loading
  if (intercept) G <- G + 1
  if (nullA) {
    A <- t(X) %*% X / nrow(X)
    A <- A[G, G, drop = F]
  }
  loading <- rep(0, p)
  loading[G] <- A %*% beta.init[G]
  loading.norm <- sqrt(sum(loading^2))

  ## Correction Direction
  direction <- compute_direction(loading, X.filter, weight.filter, deriv.filter, mu, verbose)

  ## Bias Correction
  est.plugin <- as.numeric(t(beta.init[G]) %*% A %*% beta.init[G])
  correction <- 2 * mean((weight * (y - pred) * X) %*% direction)
  est.debias <- max(as.numeric(est.plugin + correction * loading.norm), 0)

  ## Compute SE and Construct CI
  V.base <- 4 * sum(((sqrt(weight^2 * cond_var) * X) %*% direction)^2) / n^2 * loading.norm^2
  if ((n > 0.9 * p) & (model == "linear")) V.base <- V.base else V.base <- rescale^2 * V.base
  if (nullA) {
    V.A <- sum((as.vector((X[, G, drop = F] %*% beta.init[G])^2) -
      as.numeric(t(beta.init[G]) %*% A %*% beta.init[G]))^2) / n^2
  } else {
    V.A <- 0
  }
  if (model == "linear") V.add <- tau / n else V.add <- tau * max(1 / n, (sparsity * log(p) / n)^2)
  V <- (V.base + V.A + V.add)
  se <- sqrt(V)

  obj <- list(
    est.plugin = est.plugin,
    est.debias = est.debias,
    se = se,
    tau = tau
  )

  class(obj) <- "QF"
  obj
}
