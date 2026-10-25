.wname <- function(j, k, nm) if (j <= k) nm[j] else sprintf("%s(t-%d)", nm[(j - 1) %% k + 1], (j - 1) %/% k)

#' Estimate the structural VAR implied by a set of parents
#'
#' Least-squares estimate of \eqn{X_t = B_0 X_t + A_1 X_{t-1} + \dots + A_p X_{t-p} + \epsilon_t}, equation by equation, where the
#' regressors of equation i are its parents in the window (contemporaneous parents form an acyclic \eqn{B_0}).
#' @param X n by k matrix; the columns of the window are centred internally, as in the search.
#' @param p lag order.
#' @param parents list of parent sets in window indices (as returned by \code{dag_exact}).
#' @param names optional series names.
#' @return a list with \code{B0} (k by k; entry [i, j] is the coefficient of \eqn{X_j(t)} in equation i), \code{A} (list of p k by k matrices),
#'   standard errors \code{B0_se} and \code{A_se}, innovation variances \code{sig2} (maximum likelihood), \code{A0} = I - B0 (so that
#'   \eqn{A_0 X_t = \sum A_l X_{t-l} + \epsilon_t}), \code{n}, \code{names}.
#' @export
svar_estimate <- function(X, p, parents, names = NULL) {
  X <- as.matrix(X); k <- ncol(X); G <- scale(window_mat(X, p), scale = FALSE); n <- nrow(G)   # window columns centred, as in the search covariance
  nm <- if (is.null(names)) paste0("X", seq_len(k)) else names
  B0 <- matrix(0, k, k, dimnames = list(nm, nm)); B0se <- B0; A <- A_se <- vector("list", p)
  for (l in seq_len(p)) A[[l]] <- A_se[[l]] <- matrix(0, k, k, dimnames = list(nm, nm))
  sig2 <- numeric(k)
  for (i in seq_len(k)) {
    P <- sort(parents[[i]]); y <- G[, i]
    if (!length(P)) { sig2[i] <- mean(y^2); next }
    Z <- G[, P, drop = FALSE]; ZtZi <- solve(crossprod(Z)); b <- ZtZi %*% crossprod(Z, y); r <- y - Z %*% b
    sig2[i] <- mean(r^2); se <- sqrt(diag(ZtZi) * sum(r^2) / (n - length(P) - 1))   # one degree of freedom for the mean
    for (q in seq_along(P)) { j <- P[q]; col <- (j - 1) %% k + 1; lg <- (j - 1) %/% k
      if (lg == 0) { B0[i, col] <- b[q]; B0se[i, col] <- se[q] } else { A[[lg]][i, col] <- b[q]; A_se[[lg]][i, col] <- se[q] } }
  }
  list(B0 = B0, A = A, B0_se = B0se, A_se = A_se, sig2 = sig2, A0 = diag(k) - B0, n = n, names = nm, p = p)
}

#' Edge table of an estimated structural VAR
#' @param fit a \code{gmsvar} object.
#' @return data frame with columns \code{from}, \code{to}, \code{lag}, \code{coef}, \code{se}, \code{t} and \code{status}
#'   (\code{"compelled"} or \code{"reversible"} for contemporaneous edges, i.e. whether the orientation is common to the whole
#'   Markov equivalence class given the time ordering; \code{"lagged"} for lagged links, whose direction is known).
#' @export
svar_edges <- function(fit) {
  k <- fit$k; nm <- fit$names; s <- fit$svar; rows <- list()
  for (i in seq_len(k)) for (j in sort(fit$parents[[i]])) {
    lg <- (j - 1) %/% k; col <- (j - 1) %% k + 1
    cf <- if (lg == 0) s$B0[i, col] else s$A[[lg]][i, col]; se <- if (lg == 0) s$B0_se[i, col] else s$A_se[[lg]][i, col]
    st <- if (lg > 0) "lagged" else if (fit$compelled[col, i]) "compelled" else "reversible"
    rows[[length(rows) + 1]] <- data.frame(from = nm[col], to = nm[i], lag = lg, coef = cf, se = se, t = cf / se, status = st)
  }
  if (!length(rows)) return(data.frame(from = character(), to = character(), lag = integer(), coef = numeric(), se = numeric(), t = numeric(), status = character()))
  do.call(rbind, rows)
}

#' @export
coef.gmsvar <- function(object, ...) object$svar[c("B0", "A", "sig2")]

#' @export
print.gmsvar <- function(x, ...) {
  cat(sprintf("gmsvar: k = %d series, lag order p = %d, initial CIG = %s, search = %s\n", x$k, x$p, x$screen, x$search))
  for (i in seq_len(x$k)) cat(sprintf("  %s(t) <- %s\n", x$names[i], if (length(x$parents[[i]])) paste(vapply(x$parents[[i]], .wname, "", k = x$k, nm = x$names), collapse = ", ") else "(none)"))
  invisible(x)
}

#' @export
summary.gmsvar <- function(object, ...) {
  e <- svar_edges(object); cc <- e[e$lag == 0, ]
  structure(list(fit = object, edges = e, n_contemporaneous = nrow(cc), n_compelled = sum(cc$status == "compelled"), n_lagged = sum(e$lag > 0)), class = "summary.gmsvar")
}

#' @export
print.summary.gmsvar <- function(x, digits = 3, ...) {
  f <- x$fit
  crit <- if (is.null(f$ebic_gamma) || f$ebic_gamma == 0) "a BIC" else sprintf("an extended-BIC (gamma = %g)", f$ebic_gamma)
  cat(sprintf("Structural VAR(%d) identified by %s search over contemporaneous orderings (%s)\n", f$p, crit, f$search))
  cat(sprintf("n = %d, k = %d; initial CIG: %s; score %.2f\n", f$svar$n, f$k, if (is.null(f$cig)) "none (all pairs allowed)" else f$screen, f$score))
  cat(sprintf("%d contemporaneous edges (%d with orientation common to the equivalence class), %d lagged links\n\n", x$n_contemporaneous, x$n_compelled, x$n_lagged))
  e <- x$edges; e$coef <- round(e$coef, digits); e$se <- round(e$se, digits); e$t <- round(e$t, 2); print(e, row.names = FALSE)
  cat("\nInnovation standard deviations:\n"); print(round(setNames(sqrt(f$svar$sig2), f$names), digits))
  cat("\nNote: orientations marked 'reversible' are not identified by the data; the reported DAG is one member of the class.\n")
  invisible(x)
}

#' Impulse responses of the estimated SVAR
#'
#' Responses of each series to a one-standard-deviation structural shock, from the estimated \eqn{A_0}, lag matrices and shock variances.
#' Responses to shocks whose contemporaneous links are \code{"reversible"} in \code{svar_edges} depend on the member of the equivalence class chosen.
#' @param fit a \code{gmsvar} object.
#' @param h largest horizon.
#' @return array [response, shock, horizon 0..h].
#' @export
svar_irf <- function(fit, h = 12) {
  s <- fit$svar; k <- fit$k; p <- s$p; P0 <- solve(s$A0); Phi <- lapply(s$A, function(a) P0 %*% a)
  Th <- vector("list", h + 1); Th[[1]] <- diag(k)
  for (t in seq_len(h)) { acc <- matrix(0, k, k); for (l in seq_len(min(t, p))) acc <- acc + Phi[[l]] %*% Th[[t - l + 1]]; Th[[t + 1]] <- acc }
  out <- array(0, c(k, k, h + 1), dimnames = list(response = s$names, shock = s$names, horizon = 0:h))
  for (t in 0:h) out[, , t + 1] <- Th[[t + 1]] %*% P0 %*% diag(sqrt(s$sig2), k)
  out
}

#' Plot the contemporaneous DAG
#'
#' Circular layout; arrows show orientations common to the whole equivalence class, plain lines show reversible edges.
#' Edge labels give the estimated coefficients.
#' @param x a \code{gmsvar} object.
#' @param labels logical; print coefficients on the edges.
#' @param ... passed to \code{plot}.
#' @export
plot.gmsvar <- function(x, labels = TRUE, ...) {
  k <- x$k; ang <- 2 * pi * (seq_len(k) - 1) / k + pi / 2; xy <- cbind(cos(ang), sin(ang))
  graphics::plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.25, 1.25), asp = 1, axes = FALSE, xlab = "", ylab = "", ...)
  e <- svar_edges(x); e <- e[e$lag == 0, ]
  for (r in seq_len(nrow(e))) {
    a <- match(e$from[r], x$names); b <- match(e$to[r], x$names); d <- xy[b, ] - xy[a, ]; u <- d / sqrt(sum(d^2))
    p0 <- xy[a, ] + 0.16 * u; p1 <- xy[b, ] - 0.16 * u
    if (e$status[r] == "compelled") graphics::arrows(p0[1], p0[2], p1[1], p1[2], length = 0.1, lwd = 1.5) else graphics::segments(p0[1], p0[2], p1[1], p1[2], lty = 2, lwd = 1.5)
    if (labels) { lx <- p0[1] + 0.35 * (p1[1] - p0[1]); ly <- p0[2] + 0.35 * (p1[2] - p0[2]); txt <- sprintf("%.2f", e$coef[r])
      graphics::rect(lx - graphics::strwidth(txt, cex = 0.65) / 2 - 0.01, ly - graphics::strheight(txt, cex = 0.65) / 2 - 0.015, lx + graphics::strwidth(txt, cex = 0.65) / 2 + 0.01, ly + graphics::strheight(txt, cex = 0.65) / 2 + 0.015, col = "white", border = NA)
      graphics::text(lx, ly, txt, cex = 0.65, col = "grey25") }
  }
  graphics::symbols(xy[, 1], xy[, 2], circles = rep(0.15, k), inches = FALSE, add = TRUE, bg = "white")
  graphics::text(xy[, 1], xy[, 2], x$names, cex = 0.7)
  invisible(x)
}

#' Graphviz (DOT) description of the estimated structural VAR
#' @param fit a \code{gmsvar} object.
#' @param lagged include lagged links.
#' @return a character string.
#' @export
svar_dot <- function(fit, lagged = FALSE) {
  e <- svar_edges(fit); L <- c("digraph SVAR {", "  rankdir=LR;")
  for (r in seq_len(nrow(e))) {
    if (e$lag[r] > 0 && !lagged) next
    src <- if (e$lag[r] == 0) e$from[r] else sprintf("%s_lag%d", e$from[r], e$lag[r])
    att <- if (e$status[r] == "reversible") "dir=none, style=dashed" else if (e$lag[r] > 0) "color=grey" else ""
    L <- c(L, sprintf("  \"%s\" -> \"%s\" [label=\"%.2f\"%s];", src, e$to[r], e$coef[r], if (nzchar(att)) paste0(", ", att) else ""))
  }
  paste(c(L, "}"), collapse = "\n")
}

#' Moving-block bootstrap stability of the identified graph
#'
#' Re-runs \code{gmsvar} on moving-block bootstrap resamples (lag order fixed at that of \code{fit}; a supplied \code{cig} is held fixed, a screen is recomputed)
#' and records how often each contemporaneous adjacency, each orientation and each lagged link appears.
#' @param fit a \code{gmsvar} object.
#' @param X the data used for \code{fit}.
#' @param B number of resamples.
#' @param block block length.
#' @param seed seed (the caller's random stream is preserved).
#' @return list with \code{adjacency} (k by k frequencies of a contemporaneous adjacency), \code{direction} (entry [a, b] = frequency with which a -> b is a \emph{compelled} arrow, i.e. identified by the data; the rest of the adjacency frequency is reversible),
#'   \code{lagged} (list over lags of k by k frequencies, entry [i, j] for X_j(t-l) -> X_i(t)).
#' @export
gmsvar_boot <- function(fit, X, B = 200, block = 8, seed = 1) {
  X <- as.matrix(X); n <- nrow(X); k <- fit$k; p <- fit$p
  if (exists(".Random.seed", envir = globalenv())) { old <- get(".Random.seed", envir = globalenv()); on.exit(assign(".Random.seed", old, envir = globalenv())) }
  set.seed(seed); dirm <- matrix(0, k, k); adjm <- matrix(0, k, k); adjsum <- matrix(0, k, k); lagm <- replicate(p, matrix(0, k, k), simplify = FALSE)
  for (b in seq_len(B)) {
    st <- sample.int(n - block + 1, ceiling(n / block), replace = TRUE); Xb <- X[unlist(lapply(st, function(s) s:(s + block - 1)))[seq_len(n)], , drop = FALSE]
    fb <- gmsvar(Xb, p = p, cig = if (fit$cig_user) fit$cig else NULL, screen = if (fit$cig_user) "none" else fit$screen, alpha = fit$alpha, search = fit$search, ebic_gamma = if (is.null(fit$ebic_gamma)) 0 else fit$ebic_gamma, restarts = if (is.null(fit$restarts)) 3 else fit$restarts)
    for (i in seq_len(k)) for (j in fb$parents[[i]]) { lg <- (j - 1) %/% k; col <- (j - 1) %% k + 1
      if (lg == 0) adjm[col, i] <- adjm[i, col] <- 1 else lagm[[lg]][i, col] <- lagm[[lg]][i, col] + 1 }
    adjsum <- adjsum + adjm; adjm[] <- 0; dirm <- dirm + fb$compelled
  }
  dimnames(dirm) <- list(fit$names, fit$names); dimnames(adjsum) <- list(fit$names, fit$names)
  list(adjacency = adjsum / B, direction = dirm / B, lagged = lapply(lagm, function(a) { dimnames(a) <- list(fit$names, fit$names); a / B }), B = B, block = block)
}
