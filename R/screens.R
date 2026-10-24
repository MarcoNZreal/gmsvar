# screens.R -- better penalised screens for the CIG step.
#  AGL: adaptive graphical lasso, entry-wise penalty lambda / |r_ij| (r = sample partial correlation), lambda by EBIC (gamma = 0.5)
#  NW : node-wise lasso (Meinshausen-Buhlmann) for each CURRENT variable on all other window variables, lambda per node by BIC, OR rule
#  ANW: adaptive node-wise lasso, weights 1/|OLS coefficient|, lambda per node by BIC, OR rule
# All work on the window correlation matrix.  Source glasso.R first (uses .lasso_cd, glasso_path).

screen_path_nodewise <- function(S, N, k, cs = 10^seq(0, -3, length.out = 30), adaptive = FALSE) {
  m <- nrow(S); out <- list(); bic_sel <- matrix(FALSE, m, m)
  nz <- array(FALSE, c(m, m, length(cs)))          # nz[i, j, l]: coefficient of j in node i at scale c_l
  for (i in seq_len(k)) {
    idx <- setdiff(seq_len(m), i); S11 <- S[idx, idx, drop = FALSE]; s12 <- S[idx, i]
    if (adaptive) { b0 <- tryCatch(solve(S11, s12), error = function(e) rep(1, length(idx))); d <- pmax(abs(b0), 1e-8) }
    else d <- rep(1, length(idx))
    S1 <- S11 * outer(d, d); s1 <- s12 * d              # scaled problem: beta = d * beta*
    lmax <- max(abs(s1)); beta <- rep(0, length(idx)); bics <- numeric(length(cs))
    for (l in seq_along(cs)) {
      beta <- .lasso_cd_R(S1, s1, rep(cs[l] * lmax, length(idx)), beta)
      b <- beta * d; rss <- 1 - 2 * sum(b * s12) + sum(b * as.vector(S11 %*% b))
      df <- sum(b != 0); bics[l] <- N * log(max(rss, 1e-12)) + df * log(N)
      nz[i, idx, l] <- b != 0
    }
    lb <- which.min(bics); bic_sel[i, ] <- nz[i, , lb]
  }
  list(nz = nz, cs = cs, sel = bic_sel | t(bic_sel))
}

# symmetric OR-rule edge set per scale and BIC-selected edge set
screen_nodewise <- function(X, lags, k, adaptive = FALSE) {
  G <- window_mat(X, lags); S <- cor(G); N <- nrow(G)
  screen_path_nodewise(S, N, k, adaptive = adaptive)
}

screen_agl <- function(X, lags, k, ngrid = 30, gamma_idx = 2) {
  G <- window_mat(X, lags); S <- cor(G); N <- nrow(G); m <- ncol(G)
  W <- -cov2cor(solve(S)); diag(W) <- 1; Wt <- 1 / pmax(abs(W), 1e-6)       # penalty weights 1/|partial correlation|
  Wt <- Wt / mean(Wt[row(Wt) != col(Wt)])
  lam0 <- max(abs(S[row(S) != col(S)]) / Wt[row(Wt) != col(Wt)])
  lams <- lam0 * 10^seq(0, -3, length.out = ngrid)
  path <- glasso_path(S, N, lams, pen = Wt, diag0 = TRUE)
  l <- which.min(sapply(path, function(z) z$ebic[gamma_idx])); E <- path[[l]]$edges
  list(sel = E | t(E), path = path)
}

# Dispatcher used by the simulations: returns a symmetric m x m logical candidate graph
screen_graph <- function(X, lags, k, method) {
  switch(method,
    "AGL" = screen_agl(X, lags, k)$sel,
    "NW"  = screen_nodewise(X, lags, k, FALSE)$sel,
    "ANW" = screen_nodewise(X, lags, k, TRUE)$sel,
    stop("unknown screen ", method))
}
