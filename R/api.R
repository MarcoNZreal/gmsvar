#' gmsvar: graphical modelling and DAG recovery for structural VARs
#'
#' Extended-BIC search over contemporaneous orderings (over all orderings for small k, local search otherwise) from the data, the equivalence class under the
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
#' @param parents expected number of contemporaneous parents of each variable after the first in a random ordering: the number is
#'   drawn from a Binomial(3, parents / 3) distribution and capped by the number of predecessors.
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
cig_stats <- function(X, lags) {
  X <- as.matrix(X); m <- ncol(X) * (lags + 1); N <- nrow(X) - lags
  if (m >= N - 1) stop(sprintf("partial correlations of the window need more observations than window variables (N = %d, m = %d)", N, m))
  cigts_stats(X, lags)
}

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

#' Full search over orderings for the contemporaneous DAG
#'
#' Dynamic programme over all orderings of the current variables (BIC or extended BIC); the cost doubles with each added current variable (about 1 s at k = 10 and 8 s at k = 12 on one core).
#' The search is exhaustive over orderings, but the parent set of each variable is chosen by backward elimination, which can miss the best subset, so the result
#' need not be the exact minimiser of the score (in the accompanying simulations the two gave the same graph in most fits and similar accuracy).
#' @param X n by k matrix.
#' @param lags lag order.
#' @param E optional symmetric candidate graph (CIG) restricting the parents.
#' @param ebic_gamma extended-BIC parameter (0 = BIC); penalty per parent \code{log(N) + 2 * ebic_gamma * log(m - 1)}.
#' @return list with \code{parents} (parents of each current variable, window indices) and \code{score}.
#' @export
dag_exact <- function(X, lags, E = NULL, ebic_gamma = 0) learn_dag(as.matrix(X), lags, ncol(X), E, ebic_gamma = ebic_gamma)

#' Heuristic BIC search over orderings (larger k)
#'
#' Local search over orderings from a greedy start and several further starts (random orders, or, with \code{restart = "perturb"}, perturbations of the incumbent). Node problems are memoised on the sorted predecessor set. It is a heuristic: it returns a local
#' optimum of the score and the consistency result for the penalised-likelihood minimiser does not apply to it.
#' @inheritParams dag_exact
#' @param restarts number of further starts.
#' @param restart \code{"random"} (default) or \code{"perturb"} (iterated local search from the incumbent).
#' @param seed seed (the caller's random stream is preserved).
#' @export
dag_heuristic <- function(X, lags, E = NULL, restarts = 3, seed = 1, restart = "random", ebic_gamma = 0) learn_dag_heur(as.matrix(X), lags, ncol(X), E, restarts, seed, restart = restart, ebic_gamma = ebic_gamma)

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
#' restrict the DAG search. In the accompanying simulations no screen improved on the unrestricted extended-BIC search (the default of \code{gmsvar}); a CIG is useful
#' to impose prior knowledge, to inspect the data and to save computation.
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

#' Two-stage pipeline: from the data to the estimated SVAR
#'
#' Starting from the data alone, chooses the lag order (unless \code{p} is given), searches for the contemporaneous DAG over orderings
#' by the extended BIC (dynamic programme over all orderings for small k, local search otherwise; parent sets by backward elimination), derives the arrows that are common to the Markov equivalence
#' class under the time ordering, and estimates the resulting structural VAR by least squares equation by equation. An initial
#' conditional independence graph (CIG) of the window is an optional restriction (\code{cig}, e.g. from \code{cig_estimate}).
#'
#' The method has two stages. Stage 1 (optional) screens the window for a sparse conditional independence graph that restricts the candidate parents;
#' stage 2 searches over orderings with the extended BIC (\code{ebic_gamma = 0.5}). By default stage 1 is skipped (\code{screen = "none"}): in the
#' accompanying simulations the extended BIC without a screen was as accurate as, or more accurate than, the BIC with or without the adaptive
#' node-wise lasso screen at every size studied (k = 3 to 50); the BIC overfitted from k = 15. For very large systems in which the window
#' dimension k(p+1) approaches the number of observations, a liberal lasso screen (\code{screen = "NW"}) saves much computation, and when it
#' exceeds the number of observations only a lasso screen can be used. A CIG can also be supplied through \code{cig} (your own, or from \code{cig_estimate});
#' \code{screen} is a shortcut for \code{cig = cig_estimate(X, p, screen, alpha)}. A supplied \code{cig} overrides \code{screen}.
#' @param X n by k matrix (or data frame) of series; column names are used as labels.
#' @param p lag order, or \code{NULL} to choose it with \code{ic} (or from the dimension of \code{cig}).
#' @param pmax largest order considered.
#' @param ic \code{"AIC"} (default), \code{"BIC"} or \code{"HQ"}. In the accompanying simulations BIC and HQ underfitted when higher-order lags were weak, which hurt the graph; AIC overfitted slightly at a small cost.
#' @param cig optional initial CIG: a logical (or 0/1) matrix of dimension k(p+1) by k(p+1) over the window (X_t, X_{t-1}, ..., X_{t-p}),
#'   or k by k, in which case it constrains the contemporaneous pairs only and every lagged pair is allowed. Only pairs involving a current variable matter.
#' @param screen stage 1 of the two-stage method, a shortcut for \code{cig_estimate}: \code{"none"} (default), \code{"auto"} (the rule of version 0.3.0: none if k <= 20, otherwise \code{"ANW"}), a rule of \code{cig_select} (e.g. \code{"BIC"}, \code{"GMTS-ua"}, \code{"SIN"}) or a penalised screen
#'   (\code{"AGL"}, \code{"ANW"}, \code{"NW"}). Liberal screens are recommended; conservative ones can cost accuracy. Ignored if \code{cig} is given.
#' @param alpha level for the test-based screens.
#' @param search \code{"auto"} (\code{dag_exact()}, the full search over orderings, if k <= 10, else heuristic; the full search takes about 1 s at k = 10 with two lags and 8 s with six), \code{"exact"} or \code{"heuristic"}.
#' @param restarts number of further starts of the heuristic search (0 = greedy start and local search only).
#' @param ebic_gamma extended-BIC parameter: the penalty per parent is \code{log(N) + 2 * ebic_gamma * log(m - 1)}, with m the window dimension; 0 gives the BIC. The default 0.5 was about as accurate as the BIC for k <= 7 and much more accurate for k >= 20 in the accompanying simulations.
#' @return an object of class \code{gmsvar}: \code{p}, \code{order} (lag-order criteria, if used), \code{cig} (the initial graph used, or \code{NULL}),
#'   \code{parents}, \code{score}, \code{svar} (estimated SVAR, see \code{svar_estimate}), \code{compelled} (k by k logical: contemporaneous
#'   arrows common to every DAG in the equivalence class), \code{names}.
#' @examples
#' set.seed(1); X <- simulate_svar(svar_example(), 600)
#' fit <- gmsvar(X, pmax = 4)
#' fit$p; fit$parents
#' summary(fit)
#' @export
gmsvar <- function(X, p = NULL, pmax = 6, ic = "AIC", cig = NULL, screen = "none", alpha = 0.05, search = "auto", restarts = 3, ebic_gamma = 0.5) {
  ic <- match.arg(ic, c("AIC", "BIC", "HQ")); search <- match.arg(search, c("auto", "exact", "heuristic"))
  nm <- if (!is.null(colnames(X))) colnames(X) else paste0("X", seq_len(ncol(X)))
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric")
  if (anyNA(X)) stop("X contains missing values; remove or impute them first")
  sdv <- apply(X, 2, stats::sd); if (any(!is.finite(sdv)) || any(sdv <= 1e-10 * max(sdv))) stop(sprintf("X has constant column(s): %s; remove them first", paste(which(!is.finite(sdv) | sdv <= 1e-10 * max(sdv)), collapse = ", ")))
  if (qr(scale(X, scale = FALSE))$rank < ncol(X)) stop("the columns of X are linearly dependent (e.g. an aggregate and all its components); remove the redundant series first")
  X <- scale(X, scale = FALSE); k <- ncol(X); ord <- NULL
  if (identical(screen, "auto")) screen <- if (k <= 20) "none" else "ANW"   # legacy rule of version 0.3.0
  if (!is.null(cig)) {
    cig <- as.matrix(cig); d <- nrow(cig)
    if (nrow(cig) != ncol(cig)) stop("cig must be a square matrix")
    if (d %% k != 0) stop("dimension of cig must be a multiple of the number of series")
    if (is.null(p)) p <- d %/% k - 1
    if (d != k && d != k * (p + 1)) stop("cig must be k by k or k(p+1) by k(p+1)")
  } else if (is.null(p)) {
    ord <- var_order(X, pmax); p <- ord[[ic]]; if (p < 1) stop("selected order is 0; supply p")
  }
  E <- NULL; m <- k * (p + 1); N <- nrow(X) - p
  test_screen <- is.null(cig) && !(screen %in% c("none", "AGL", "ANW", "NW"))
  if (test_screen && m >= N - 1)
    stop(sprintf("partial-correlation screens need more observations than window variables (N = %d, m = %d); use screen = \"NW\" or \"ANW\"", N, m))
  if (is.null(cig) && screen == "none" && m - 1 >= N - 1)
    warning(sprintf("the window has m = %d variables but only N = %d observations: the unscreened search is not reliable; consider a lasso screen (screen = \"NW\")", m, N))
  if (!is.null(cig)) {
    if (nrow(cig) == k && m > k) { E <- matrix(TRUE, m, m); E[seq_len(k), seq_len(k)] <- cig != 0 } else E <- cig != 0
    E <- E | t(E)
  } else if (screen %in% c("AGL", "ANW", "NW")) E <- cig_screen(X, p, screen)
  else if (screen != "none") E <- cig_select(cig_stats(X, p), k, p, screen, alpha)
  if (search == "auto") search <- if (k <= 10) "exact" else "heuristic"
  fit <- if (search == "exact") dag_exact(X, p, E, ebic_gamma = ebic_gamma) else dag_heuristic(X, p, E, restarts = restarts, ebic_gamma = ebic_gamma)
  out <- list(p = p, order = ord, cig = E, parents = fit$parents, score = fit$score, k = k, names = nm,
              screen = if (is.null(cig)) screen else "user cig", search = search, cig_user = !is.null(cig), alpha = alpha, ebic_gamma = ebic_gamma, restarts = restarts)
  out$svar <- svar_estimate(X, p, fit$parents, nm)
  out$compelled <- cpdag_cc(fit$parents, k, k * (p + 1))
  structure(out, class = "gmsvar")
}
