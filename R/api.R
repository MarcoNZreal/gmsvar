#' gmsvar: graphical modelling and DAG recovery for structural VARs
#'
#' BIC search over contemporaneous orderings (exact for small k, heuristic otherwise) from the data, the equivalence class under the
#' time ordering, and the estimated SVAR; window partial-correlation graphs (CIGs), edge-selection rules and penalised screens are optional tools.
#' @keywords internal
#' @importFrom stats cor cov cov2cor pnorm pt rbinom rnorm rt runif coef setNames
#' @importFrom graphics plot arrows segments symbols text
#' @importFrom utils combn
"_PACKAGE"

#' Example SVAR(2) in three variables
#' @return a list with \code{A0}, \code{A} (list of lag matrices) and \code{sig2}.
#' @export
svar_example <- function() model1()

#' Random sparse stable SVAR
#'
#' Draws a sparse acyclic \eqn{A_0}, lag matrices with autoregressive diagonals and sparse off-diagonal entries,
#' and innovation variances, rescaling the lag matrices until the largest companion root is below \code{max_root}.
#' The call is reproducible from \code{seed} and leaves the caller's random-number stream unchanged.
#' @param k number of variables.
#' @param p lag order.
#' @param seed seed.
#' @param max_root bound on the largest companion root.
#' @param parents expected number of contemporaneous parents per variable (before capping by the ordering).
#' @param off1 number of non-zero off-diagonal entries of the lag-1 matrix.
#' @param off2 number of non-zero entries of each lag matrix of order 2 and higher.
#' @return an SVAR list (see \code{svar_example}).
#' @export
random_svar <- function(k = 15, p = 2, seed = 15015, max_root = 0.9, parents = 1.3, off1 = round(1.5 * k), off2 = round(0.93 * k))
  random_svar_impl(k, p, seed, parents = parents, off1 = off1, off2 = off2, max_root = max_root)

#' Simulate a Gaussian (or Student t) SVAR
#' @param m an SVAR list.
#' @param n sample size.
#' @param burn burn-in length.
#' @param errors \code{"gaussian"} or \code{"t"}.
#' @param df degrees of freedom for t errors.
#' @return an n by k matrix.
#' @export
simulate_svar <- function(m, n, burn = 500, errors = "gaussian", df = 5) simulate_svar_impl(m, n, burn, errors, df)

#' Window matrix
#' @param X n by k matrix.
#' @param lags lag order.
#' @return matrix with columns X_t, X_{t-1}, ..., X_{t-lags} (current variables first).
#' @export
window_matrix <- function(X, lags) window_mat(as.matrix(X), lags)

#' Partial correlations and tests in the window
#' @param X n by k matrix.
#' @param lags lag order.
#' @return list with \code{R} (partial correlations), \code{t}, \code{pt} (t-test p-values), \code{pz} (Fisher-z p-values), \code{nv}, \code{m}.
#' @export
cig_stats <- function(X, lags) cigts_stats(as.matrix(X), lags)

#' Select CIG edges from partial-correlation tests
#'
#' All rules rank the pairs identically (by the absolute partial correlation); they differ only in the cutoff.
#' @param stats output of \code{cig_stats}.
#' @param k number of series.
#' @param lags lag order.
#' @param method one of \code{"GMTS-ua"} (t-test at level alpha), \code{"GMTS-Bonf"}, \code{"GMTS-Sidak"}, \code{"SIN"}
#'   (Sidak-Holm step-down on Fisher-z p-values), \code{"BH"}, \code{"BIC"} (cutoff \eqn{|t| > \sqrt{\log n}}).
#' @param alpha significance level.
#' @return symmetric logical m by m matrix of selected edges.
#' @export
cig_select <- function(stats, k, lags, method = "BIC", alpha = 0.05) { E <- select_edges(stats, k, lags, method, alpha); E | t(E) }

#' Penalised screens for the CIG
#' @param X n by k matrix.
#' @param lags lag order.
#' @param method \code{"AGL"} (adaptive graphical lasso, EBIC), \code{"ANW"} (adaptive node-wise lasso) or \code{"NW"} (node-wise lasso).
#' @return symmetric logical m by m matrix.
#' @export
cig_screen <- function(X, lags, method = c("ANW", "AGL", "NW")) screen_graph(as.matrix(X), lags, ncol(X), match.arg(method))

#' Exact BIC search for the contemporaneous DAG
#'
#' Dynamic programme over subsets of current variables; feasible for about k <= 8.
#' @param X n by k matrix.
#' @param lags lag order.
#' @param E optional symmetric candidate graph (CIG) restricting the parents.
#' @return list with \code{parents} (parents of each current variable, window indices) and \code{score}.
#' @export
dag_exact <- function(X, lags, E = NULL) learn_dag(as.matrix(X), lags, ncol(X), E)

#' Heuristic BIC search over orderings (larger k)
#'
#' Local search over orderings from a greedy and several random starts. It is a heuristic: it returns a local
#' optimum of the BIC and the consistency result for the exact search does not apply to it.
#' @inheritParams dag_exact
#' @param restarts number of random starts.
#' @param seed seed (the caller's random stream is preserved).
#' @export
dag_heuristic <- function(X, lags, E = NULL, restarts = 3, seed = 1) learn_dag_heur(as.matrix(X), lags, ncol(X), E, restarts, seed)

#' Parent sets of a model and DAG error metrics
#' @param m an SVAR list.
#' @return a list of true parent sets in window indices.
#' @export
true_parent_sets <- function(m) true_parents(m)

#' Compare an estimated DAG with the truth
#' @param est,truth lists of parent sets.
#' @param k number of series.
#' @param p lag order.
#' @return named vector with skeleton counts (contemporaneous and lagged) and indicators of exact skeleton and equivalence-class recovery.
#' @export
dag_compare <- function(est, truth, k, p) dag_metrics(est, truth, k, k * (p + 1))

#' Estimate an initial CIG (optional step)
#'
#' Convenience wrapper that returns a logical k(p+1) by k(p+1) conditional independence graph of the window, by a partial-correlation
#' test rule (\code{cig_select}) or a penalised screen (\code{cig_screen}). Pass the result to \code{gmsvar(X, p, cig = ...)} to
#' restrict the DAG search. The accompanying simulations found that the unrestricted search was the most accurate; a CIG is useful
#' to impose prior knowledge, to inspect the data, or (liberal screens only) to speed up the search when k is large.
#' @param X n by k matrix.
#' @param p lag order.
#' @param method a rule of \code{cig_select} (\code{"GMTS-ua"}, \code{"GMTS-Bonf"}, \code{"GMTS-Sidak"}, \code{"SIN"}, \code{"BH"}, \code{"BIC"})
#'   or a penalised screen (\code{"AGL"}, \code{"ANW"}, \code{"NW"}).
#' @param alpha level for the test-based rules.
#' @return symmetric logical matrix of dimension k(p+1), columns ordered X_t, X_{t-1}, ..., X_{t-p}.
#' @examples
#' set.seed(1); X <- simulate_svar(svar_example(), 600)
#' E <- cig_estimate(X, 2, "GMTS-ua", alpha = 0.3)
#' fit <- gmsvar(X, p = 2, cig = E)
#' @export
cig_estimate <- function(X, p, method = "GMTS-ua", alpha = 0.05) {
  X <- scale(as.matrix(X), scale = FALSE); k <- ncol(X)
  if (method %in% c("AGL", "ANW", "NW")) cig_screen(X, p, method) else cig_select(cig_stats(X, p), k, p, method, alpha)
}

#' End-to-end pipeline: from the data to the estimated SVAR
#'
#' Starting from the data alone, chooses the lag order (unless \code{p} is given), searches for the contemporaneous DAG over orderings
#' by the BIC (exact dynamic programme for small k, heuristic otherwise), derives the arrows that are common to the Markov equivalence
#' class under the time ordering, and estimates the resulting structural VAR by least squares equation by equation. An initial
#' conditional independence graph (CIG) of the window is an optional restriction (\code{cig}, e.g. from \code{cig_estimate}).
#'
#' The default uses no CIG: the search considers every pair, which was the most accurate option in the accompanying simulations
#' whenever it is feasible. A CIG can be supplied through \code{cig} (your own, or from \code{cig_estimate}); the argument \code{screen}
#' is a shortcut for \code{cig = cig_estimate(X, p, screen, alpha)}. A supplied \code{cig} overrides \code{screen}.
#' @param X n by k matrix (or data frame) of series; column names are used as labels.
#' @param p lag order, or \code{NULL} to choose it with \code{ic} (or from the dimension of \code{cig}).
#' @param pmax largest order considered.
#' @param ic \code{"AIC"} (default), \code{"BIC"} or \code{"HQ"}. In the accompanying simulations BIC and HQ underfitted when higher-order lags were weak, which hurt the graph; AIC overfitted slightly at a small cost.
#' @param cig optional initial CIG: a logical (or 0/1) matrix of dimension k(p+1) by k(p+1) over the window (X_t, X_{t-1}, ..., X_{t-p}),
#'   or k by k, in which case it constrains the contemporaneous pairs only and every lagged pair is allowed. Only pairs involving a current variable matter.
#' @param screen shortcut for \code{cig_estimate}: \code{"none"} (default), a rule of \code{cig_select} (e.g. \code{"BIC"}, \code{"GMTS-ua"}, \code{"SIN"}) or a penalised screen
#'   (\code{"AGL"}, \code{"ANW"}, \code{"NW"}). Liberal screens are recommended; conservative ones can cost accuracy. Ignored if \code{cig} is given.
#' @param alpha level for the test-based screens.
#' @param search \code{"auto"} (exact if k <= 7, else heuristic), \code{"exact"} or \code{"heuristic"}.
#' @return an object of class \code{gmsvar}: \code{p}, \code{order} (lag-order criteria, if used), \code{cig} (the initial graph used, or \code{NULL}),
#'   \code{parents}, \code{score}, \code{svar} (estimated SVAR, see \code{svar_estimate}), \code{compelled} (k by k logical: contemporaneous
#'   arrows common to every DAG in the equivalence class), \code{names}.
#' @examples
#' set.seed(1); X <- simulate_svar(svar_example(), 600)
#' fit <- gmsvar(X, pmax = 4)
#' fit$p; fit$parents
#' summary(fit)
#' @export
gmsvar <- function(X, p = NULL, pmax = 6, ic = "AIC", cig = NULL, screen = "none", alpha = 0.05, search = "auto") {
  nm <- if (!is.null(colnames(X))) colnames(X) else paste0("X", seq_len(ncol(X)))
  X <- scale(as.matrix(X), scale = FALSE); k <- ncol(X); ord <- NULL
  if (!is.null(cig)) {
    cig <- as.matrix(cig); d <- nrow(cig)
    if (nrow(cig) != ncol(cig)) stop("cig must be a square matrix")
    if (d %% k != 0) stop("dimension of cig must be a multiple of the number of series")
    if (is.null(p)) p <- d %/% k - 1
    if (d != k && d != k * (p + 1)) stop("cig must be k by k or k(p+1) by k(p+1)")
  } else if (is.null(p)) {
    ord <- var_order(X, pmax); p <- ord[[ic]]; if (p < 1) stop("selected order is 0; supply p")
  }
  E <- NULL; m <- k * (p + 1)
  if (!is.null(cig)) {
    if (nrow(cig) == k && m > k) { E <- matrix(TRUE, m, m); E[seq_len(k), seq_len(k)] <- cig != 0 } else E <- cig != 0
    E <- E | t(E)
  } else if (screen %in% c("AGL", "ANW", "NW")) E <- cig_screen(X, p, screen)
  else if (screen != "none") E <- cig_select(cig_stats(X, p), k, p, screen, alpha)
  if (search == "auto") search <- if (k <= 7) "exact" else "heuristic"
  fit <- if (search == "exact") dag_exact(X, p, E) else dag_heuristic(X, p, E)
  out <- list(p = p, order = ord, cig = E, parents = fit$parents, score = fit$score, k = k, names = nm,
              screen = if (is.null(cig)) screen else "user cig", search = search, cig_user = !is.null(cig), alpha = alpha)
  out$svar <- svar_estimate(X, p, fit$parents, nm)
  out$compelled <- cpdag_cc(fit$parents, k, k * (p + 1))
  structure(out, class = "gmsvar")
}
