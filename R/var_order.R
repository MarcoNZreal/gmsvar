#' Select the lag order of a VAR by information criteria
#'
#' Fits VAR(p) by least squares for p = 0, ..., pmax on a common sample (the last n - pmax observations, with
#' an intercept) and returns AIC, BIC (Schwarz) and Hannan-Quinn criteria
#' \eqn{\log\det\hat\Sigma_p + c_N k^2 p / N}, with \eqn{c_N = 2, \log N, 2\log\log N}.
#'
#' @param X numeric matrix, n observations by k series (columns).
#' @param pmax largest lag order considered (default 8).
#' @return A list with \code{table} (criterion values by p), \code{AIC}, \code{BIC}, \code{HQ} (selected orders).
#' @examples
#' set.seed(1); X <- simulate_svar(svar_example(), 400)
#' var_order(X, pmax = 5)$BIC
#' @export
var_order <- function(X, pmax = 8) {
  X <- as.matrix(X); n <- nrow(X); k <- ncol(X); N <- n - pmax
  if (N <= k * (pmax + 1)) stop("sample too short for pmax")
  Y <- X[(pmax + 1):n, , drop = FALSE]
  crit <- t(sapply(0:pmax, function(p) {
    Z <- if (p == 0) matrix(1, N, 1) else cbind(1, do.call(cbind, lapply(seq_len(p), function(l) X[(pmax + 1 - l):(n - l), , drop = FALSE])))
    R <- Y - Z %*% solve(crossprod(Z), crossprod(Z, Y))
    ld <- as.numeric(determinant(crossprod(R) / N, logarithm = TRUE)$modulus); np <- k * k * p
    c(p = p, logdet = ld, AIC = ld + 2 * np / N, BIC = ld + log(N) * np / N, HQ = ld + 2 * log(log(N)) * np / N)
  }))
  list(table = crit, AIC = unname(crit[which.min(crit[, "AIC"]), "p"]), BIC = unname(crit[which.min(crit[, "BIC"]), "p"]), HQ = unname(crit[which.min(crit[, "HQ"]), "p"]))
}
