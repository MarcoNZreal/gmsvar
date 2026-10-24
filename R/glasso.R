# glasso.R -- graphical lasso (Friedman, Hastie & Tibshirani 2008) in base R, with EBIC selection.
# Block coordinate descent; checked against the glasso package (same edges and precision matrix).

.lasso_cd_R <- function(W11, s12, lam, beta, maxit = 200, tol = 1e-7) {
  # min 0.5 b'W11 b - b's12 + lam |b|_1 by coordinate descent, warm start beta
  q <- length(s12)
  for (it in seq_len(maxit)) {
    delta <- 0
    for (i in seq_len(q)) {
      r <- s12[i] - sum(W11[i, ] * beta) + W11[i, i] * beta[i]
      bn <- sign(r) * max(abs(r) - lam[i], 0) / W11[i, i]
      delta <- max(delta, abs(bn - beta[i])); beta[i] <- bn
    }
    if (delta < tol) break
  }
  beta
}

# lam: scalar (uniform penalty) or p x p matrix of entry-wise penalties
glasso_fit <- function(S, lam, init = NULL, maxit = 100, tol = 1e-5) {
  p <- nrow(S)
  Lam <- if (length(lam) == 1) matrix(lam, p, p) else lam
  W <- if (is.null(init)) S + diag(diag(Lam), p) else init$W
  B <- if (is.null(init)) matrix(0, p - 1, p) else init$B
  for (it in seq_len(maxit)) {
    W0 <- W
    for (j in seq_len(p)) {
      B[, j] <- .lasso_cd_R(W[-j, -j, drop = FALSE], S[-j, j], Lam[-j, j], B[, j])
      w12 <- W[-j, -j, drop = FALSE] %*% B[, j]
      W[-j, j] <- w12; W[j, -j] <- w12
    }
    if (mean(abs(W - W0)) < tol * mean(abs(S[row(S) != col(S)]))) break
  }
  Theta <- matrix(0, p, p)
  for (j in seq_len(p)) {
    t22 <- 1 / (W[j, j] - sum(W[-j, j] * B[, j]))
    Theta[j, j] <- t22; Theta[-j, j] <- -B[, j] * t22
  }
  Theta <- (Theta + t(Theta)) / 2
  nz <- (abs(Theta) > 1e-10) & (row(Theta) != col(Theta))
  list(W = W, B = B, Theta = Theta, edges = nz, iter = it)
}

# Path over a decreasing lambda grid with warm starts; returns edges per lambda and EBIC values.
# pen: optional 0/1 matrix; entry (i,j) is penalised with lam * pen[i,j] (diagonal always lam)
glasso_path <- function(S, N, lams, gammas = c(0, 0.5), pen = NULL, diag0 = FALSE) {
  p <- nrow(S); out <- vector("list", length(lams)); fit <- NULL
  for (l in seq_along(lams)) {
    L <- if (is.null(pen)) lams[l] else { P <- pen * lams[l]; diag(P) <- if (diag0) 0 else lams[l]; P }
    fit <- glasso_fit(S, L, init = fit)
    ld <- as.numeric(determinant(fit$Theta, logarithm = TRUE)$modulus)
    ll <- 0.5 * N * (ld - sum(S * fit$Theta))
    E <- sum(fit$edges) / 2
    out[[l]] <- list(edges = fit$edges, ebic = vapply(gammas, function(g) -2 * ll + E * log(N) + 4 * g * E * log(p), 0))
  }
  out
}

# Candidate graph for the DAG search: glasso path, lambda chosen by (E)BIC; returns symmetric m x m logical matrix.
glasso_graph <- function(X, lags, gamma_idx = 2, variant = "U", k = ncol(X), ngrid = 30) {
  G <- window_mat(X, lags); N <- nrow(G); S <- cor(G); mm <- ncol(G)
  pen <- NULL
  if (variant == "S") { pen <- matrix(1, mm, mm); pen[(k + 1):mm, (k + 1):mm] <- 0 }
  lams <- max(abs(S[row(S) != col(S)])) * 10^seq(0, -3, length.out = ngrid)
  path <- glasso_path(S, N, lams, pen = pen)
  l <- which.min(sapply(path, function(z) z$ebic[gamma_idx]))
  E <- path[[l]]$edges; E | t(E)
}

