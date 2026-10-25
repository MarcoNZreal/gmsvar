# sim.R -- SVAR simulation, window matrix, partial-correlation statistics (port of cigts), edge-selection rules.
# Window variables are indexed 1..m: current variables 1..k, then the lag-1 block, ...
model1 <- function() {
  mat <- function(v, nr) matrix(v, nrow = nr, byrow = TRUE)
  list(name = "M1",
       A0 = mat(c(1, 0, 0,  -1, 1, 0,  1, 0, 1), 3),
       A  = list(mat(c(.9, 0, 0,  0, .6, 0,  0, 0, -.5), 3),
                 mat(c(-.6, 0, 0,  .9, 0, .4,  0, 0, 0), 3)),
       sig2 = c(.09, .04, .025))
}

random_svar_impl <- function(k = 15, p = 2, seed = 15015, parents = 1.3, off1 = round(1.5 * k), off2 = round(k * 0.93), max_root = 0.9) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed); ord <- sample(k)
  A0 <- diag(k)
  for (r in 2:k) { i <- ord[r]; np <- min(r - 1, rbinom(1, 3, parents / 3)); if (np > 0) { j <- ord[sample(r - 1, np)]; A0[i, j] <- -sample(c(-1, 1), np, TRUE) * runif(np, .3, .7) } }
  A <- vector("list", p)
  A[[1]] <- diag(runif(k, .35, .6)); offd <- setdiff(seq_len(k * k), seq(1, k * k, by = k + 1))
  idx <- sample(offd, min(off1, length(offd))); A[[1]][idx] <- sample(c(-1, 1), length(idx), TRUE) * runif(length(idx), .15, .35)
  if (p >= 2) for (l in 2:p) { Al <- matrix(0, k, k); idx <- sample(seq_len(k * k), min(off2, k * k)); Al[idx] <- sample(c(-1, 1), length(idx), TRUE) * runif(length(idx), .1, .25); A[[l]] <- Al }
  m <- list(name = sprintf("random k=%d p=%d seed=%d", k, p, seed), A0 = A0, A = A, sig2 = runif(k, .5, 1.5))
  repeat { if (max_modulus_root_inv(m) < max_root) break; m$A <- lapply(m$A, function(a) a * 0.9) }
  m
}

reduced_form <- function(m) {
  A0inv <- solve(m$A0)
  list(Phi = lapply(m$A, function(Ai) A0inv %*% Ai),
       Se = A0inv %*% diag(m$sig2) %*% t(A0inv), A0inv = A0inv)
}

companion <- function(Phi) {
  k <- nrow(Phi[[1]]); p <- length(Phi)
  F <- matrix(0, k * p, k * p)
  F[1:k, ] <- do.call(cbind, Phi)
  if (p > 1) F[(k + 1):(k * p), 1:(k * (p - 1))] <- diag(k * (p - 1))
  F
}

# largest eigenvalue modulus of the companion matrix (< 1  <=>  stationary)
max_modulus_root_inv <- function(m) max(Mod(eigen(companion(reduced_form(m)$Phi), only.values = TRUE)$values))

simulate_svar_impl <- function(m, n, burn = 500, errors = "gaussian", df = 5) {
  rf <- reduced_form(m); k <- nrow(m$A0); p <- length(m$A)
  N <- n + burn
  eps <- if (errors == "t") rt(N * k, df) * sqrt((df - 2) / df) else rnorm(N * k)
  eps <- sweep(matrix(eps, N, k), 2, sqrt(m$sig2), "*")
  e <- eps %*% t(rf$A0inv)
  Phic <- do.call(cbind, rf$Phi)
  X <- matrix(0, N, k); s <- numeric(k * p)
  for (t in seq_len(N)) {
    xn <- as.vector(Phic %*% s) + e[t, ]
    X[t, ] <- xn
    s <- c(xn, s[seq_len(k * (p - 1))])
  }
  X[(burn + 1):N, , drop = FALSE]
}

# ------------------------------------------------------------------ cigts
# window matrix: columns [X_t, X_{t-1}, ..., X_{t-lags}], rows t = lags+1..n
window_mat <- function(X, lags) {
  n <- nrow(X)
  do.call(cbind, lapply(0:lags, function(i) X[(lags + 1 - i):(n - i), , drop = FALSE]))
}

# cigts.m: V = corrcoef(G); W = inv(V); R = -W / (s s');  df = nv - (m-1)
cigts_stats <- function(X, lags) {
  G <- window_mat(X, lags)
  nv <- nrow(G); m <- ncol(G)
  W <- solve(cor(G))
  s <- sqrt(diag(W))
  R <- -W / outer(s, s); diag(R) <- 1
  df <- nv - (m - 1)
  r <- pmax(pmin(R, 1 - 1e-9), -1 + 1e-9)
  tt <- r * sqrt(df) / sqrt(1 - r^2)
  pt_ <- 2 * pt(abs(tt), df, lower.tail = FALSE)
  z <- atanh(r) * sqrt(nv - (m - 2) - 3)
  pz <- 2 * pnorm(abs(z), lower.tail = FALSE)
  list(R = R, df = df, nv = nv, m = m, pt = pt_, pz = pz, t = tt)
}

crit_value <- function(tv, df) tv / sqrt(tv^2 + df)   # eq. (3.4) / cigts.m

# cells (i<j) with at least one current variable (i<=k): excludes the trivial diagonal
scored_mask <- function(k, lags) {
  m <- k * (lags + 1); M <- matrix(FALSE, m, m)
  for (i in 1:k) if (i < m) M[i, (i + 1):m] <- TRUE
  M
}
all_pairs_mask <- function(m) upper.tri(matrix(0, m, m))

# ------------------------------------------------------------ multiplicity
holm_sidak <- function(p) {
  M <- length(p); o <- order(p)
  adj <- cummax(1 - (1 - p[o])^(M - (0:(M - 1))))
  out <- numeric(M); out[o] <- pmin(adj, 1); out
}
bh_adjust <- function(p) {
  M <- length(p); o <- order(p)
  q <- rev(cummin(rev(p[o] * M / seq_len(M))))
  out <- numeric(M); out[o] <- pmin(q, 1); out
}

# Returns an m x m logical (upper-triangular) edge matrix over all window pairs.
select_edges <- function(S, k, lags, method, alpha = 0.05) {
  m <- S$m; U <- all_pairs_mask(m)
  pt_ <- S$pt[U]; pz <- S$pz[U]; M <- sum(U)
  sel <- switch(method,
    "GMTS-ua"    = pt_ < alpha,
    "GMTS-Bonf"  = pt_ < alpha / M,
    "GMTS-BonfS" = pt_ < alpha / sum(scored_mask(k, lags)),
    "GMTS-Sidak" = pt_ < 1 - (1 - alpha)^(1 / M),
    "SIN"        = holm_sidak(pz) < alpha,
    "BH"         = bh_adjust(pt_) < alpha,
    "BIC"        = abs(S$t[U]) > sqrt(log(S$nv)),
    stop("unknown method: ", method))
  E <- matrix(FALSE, m, m); E[U] <- sel; E
}

