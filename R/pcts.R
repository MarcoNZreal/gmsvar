# pc_ts.R -- tiered PC-stable learner for SVAR window data (PCMCI+ style).
# * Gaussian partial-correlation (Fisher z) tests on the window vector; lag-lag pairs are never tested (always adjacent, as in the
#   augmented graph); lag -> current edges are oriented by the arrow of time (tiers).
# * Conditioning sets are drawn from the (stable) adjacencies of the two endpoints, ordered by association strength (largest
#   |partial correlation| first) and capped at the maxc strongest candidates and subsets of size <= maxl, as in the PC1 stage of PCMCI.
# * Orientation: v-structures from the separating sets, then Meek rules R1-R3; the PDAG is extended to a DAG (Dor-Tarsi).
# alpha_pc >= 1 disables the PC1 pruning (all lagged variables are kept as conditioning set; the default here, because PC1 pruning with lag-only
# conditioning lost many true lagged parents in our designs).
# This is a PCMCI+-style learner written from scratch for this study, not the tigramite implementation (no MCI stage).

# all subsets of size l of x (safe for length(x) == 1, unlike combn)
subsets <- function(x, l) if (length(x) == l) list(x) else combn(x, l, simplify = FALSE)

pc_ts_impl <- function(X, lags, k, alpha = 0.01, alpha_pc = 1.01, maxl = 3, maxc = 6, maxpc = 3, use_shift = FALSE) {
  G <- window_mat(X, lags); N <- nrow(G); m <- ncol(G); S <- cor(G); p <- lags
  Om <- solve(S); Rp <- abs(-Om / sqrt(outer(diag(Om), diag(Om))))
  adj <- matrix(TRUE, m, m); diag(adj) <- FALSE
  sep <- new.env(hash = TRUE); key <- function(a, b) paste(min(a, b), max(a, b))
  pval <- function(a, b, cond) {
    cond <- unique(setdiff(cond, c(a, b))); idx <- c(a, b, cond)
    Oi <- tryCatch(solve(S[idx, idx, drop = FALSE]), error = function(e) NULL)
    if (is.null(Oi)) return(1)
    r <- -Oi[1, 2] / sqrt(Oi[1, 1] * Oi[2, 2]); r <- max(min(r, 1 - 1e-12), -1 + 1e-12)
    df <- N - length(cond) - 3; if (df < 1) return(1)
    2 * pnorm(abs(atanh(r)) * sqrt(df), lower.tail = FALSE)
  }
  lagv <- if (m > k) (k + 1):m else integer(0)
  # ---- Stage 1 (PC1): lagged parent candidates of each current variable, conditioning on lagged variables only
  Bset <- vector("list", k)
  for (i in seq_len(k)) {
    P <- lagv
    for (l in 0:maxpc) {
      if (length(P) <= l) break
      st <- P; rem <- integer(0)
      for (cnd in st) {
        oth <- setdiff(st, cnd); if (length(oth) > maxc) oth <- oth[order(Rp[i, oth], decreasing = TRUE)][seq_len(maxc)]
        subs <- if (l == 0) list(integer(0)) else subsets(oth, l)
        for (cs in subs) if (pval(i, cnd, cs) > alpha_pc) { rem <- c(rem, cnd); break }
      }
      P <- setdiff(P, rem)
    }
    Bset[[i]] <- P[order(Rp[i, P], decreasing = TRUE)]
  }
  shifted <- function(cl) {                 # lagged parents of the lag vertex cl, shifted into the window (truncated at lag p)
    tau <- (cl - 1) %/% k; j <- (cl - 1) %% k + 1; out <- integer(0)
    for (q in Bset[[j]]) { tq <- (q - 1) %/% k; l2 <- (q - 1) %% k + 1; if (tau + tq <= p) out <- c(out, (tau + tq) * k + l2) }
    out
  }
  # lag-lag pairs and non-candidate lagged links
  for (i in seq_len(k)) for (cl in lagv) if (!(cl %in% Bset[[i]])) { adj[i, cl] <- adj[cl, i] <- FALSE; assign(key(i, cl), Bset[[i]], envir = sep) }
  # ---- Stage 2: MCI-type tests for lag -> current links and PC-stable tests for contemporaneous links
  for (l in 0:maxl) {
    snap <- adj
    for (i in seq_len(k)) {
      for (cl in lagv) if (adj[i, cl]) {
        cn <- setdiff(which(snap[i, seq_len(k)]), i)
        if (length(cn) > maxc) cn <- cn[order(Rp[i, cn], decreasing = TRUE)][seq_len(maxc)]
        if (length(cn) < l) next
        base <- if (use_shift) c(setdiff(Bset[[i]], cl), shifted(cl)) else setdiff(Bset[[i]], cl)
        subs <- if (l == 0) list(integer(0)) else subsets(cn, l)
        for (cs in subs) if (pval(i, cl, c(base, cs)) > alpha) { adj[i, cl] <- adj[cl, i] <- FALSE; assign(key(i, cl), c(base, cs), envir = sep); break }
      }
    }
    for (i in seq_len(k - 1)) for (j in (i + 1):k) if (adj[i, j]) {
      base <- union(Bset[[i]], Bset[[j]]); done <- FALSE
      for (ends in list(c(i, j), c(j, i))) {
        x <- ends[1]; y <- ends[2]
        cn <- setdiff(which(snap[x, seq_len(k)]), c(x, y))
        if (length(cn) > maxc) cn <- cn[order(Rp[x, cn], decreasing = TRUE)][seq_len(maxc)]
        if (length(cn) < l) next
        subs <- if (l == 0) list(integer(0)) else subsets(cn, l)
        for (cs in subs) if (pval(i, j, c(base, cs)) > alpha) { adj[i, j] <- adj[j, i] <- FALSE; assign(key(i, j), c(base, cs), envir = sep); done <- TRUE; break }
        if (done) break
      }
    }
  }
  # ---- orientation
  D <- matrix(FALSE, m, m)                                   # D[a, b]: a -> b
  for (a in seq_len(m)) for (b in seq_len(m)) if (adj[a, b]) {
    if (a > k && b <= k) D[a, b] <- TRUE                     # lag -> current
    if (a > k && b > k && a < b) D[a, b] <- TRUE             # lag-lag fixed by index
  }
  for (cc in seq_len(k)) {
    nb <- which(adj[cc, ])
    if (length(nb) < 2) next
    for (ii in seq_len(length(nb) - 1)) for (jj in (ii + 1):length(nb)) {
      a <- nb[ii]; b <- nb[jj]
      if (adj[a, b]) next
      ss <- get0(key(a, b), envir = sep, inherits = FALSE)
      if (!is.null(ss) && !(cc %in% ss)) {
        if (!D[cc, a]) D[a, cc] <- TRUE
        if (!D[cc, b]) D[b, cc] <- TRUE
      }
    }
  }
  und <- function() adj & !D & !t(D)
  repeat {                                                    # Meek R1-R3
    ch <- FALSE; U <- und()
    for (b in seq_len(k)) for (cc in seq_len(k)) if (U[b, cc]) {
      if (any(D[, b] & !adj[, cc] & seq_len(m) != cc)) { D[b, cc] <- TRUE; ch <- TRUE; U <- und(); next }                  # R1
      if (any(D[b, ] & D[, cc])) { D[b, cc] <- TRUE; ch <- TRUE; U <- und(); next }                                       # R2
      cand <- which(U[b, ] & D[, cc] & seq_len(m) != cc)                                                                  # R3
      if (length(cand) >= 2) { pr <- combn(cand, 2); if (any(apply(pr, 2, function(z) !adj[z[1], z[2]]))) { D[b, cc] <- TRUE; ch <- TRUE; U <- und() } }
    }
    if (!ch) break
  }
  # ---- DAG extension (Dor-Tarsi); fallback orients by index
  Dd <- D; Ud <- und(); alive <- rep(TRUE, m); progress <- TRUE
  while (any(alive) && progress) {
    progress <- FALSE
    for (x in which(alive)) {
      if (any(Dd[x, alive])) next
      nb_u <- which(Ud[x, ] & alive); ad_x <- which((adj[x, ]) & alive)
      ok <- TRUE
      for (y in nb_u) { need <- setdiff(ad_x, y); if (length(need) && !all(adj[y, need])) { ok <- FALSE; break } }
      if (!ok) next
      for (y in nb_u) { Dd[y, x] <- TRUE; Ud[x, y] <- Ud[y, x] <- FALSE }
      alive[x] <- FALSE; progress <- TRUE; break
    }
  }
  Ud <- adj & !Dd & !t(Dd)
  for (a in seq_len(m)) for (b in seq_len(m)) if (a < b && Ud[a, b]) Dd[b, a] <- TRUE   # fallback
  parents <- lapply(seq_len(k), function(i) sort(which(Dd[, i])))
  list(parents = parents, adj = adj, D = D, alpha = alpha, sep = as.list(sep))
}

#' Tiered PC-stable learner for the contemporaneous DAG (PCMCI+-style, written for this package)
#'
#' Constraint-based comparator to the BIC search: Gaussian partial-correlation (Fisher z) tests on the window, lag-lag pairs always
#' adjacent, lag-to-current links oriented by the arrow of time, v-structures from separating sets, Meek rules and a Dor-Tarsi
#' extension to a DAG. The PC1 pruning stage is switched off by default. This is not the \pkg{tigramite} implementation
#' (there is no MCI stage), and in the accompanying simulations it was never more accurate than the BIC search.
#' @inheritParams dag_exact
#' @param alpha significance level of the tests.
#' @param maxl,maxc maximum conditioning-set size and number of candidate conditioning variables.
#' @return list with \code{parents} (as for \code{dag_exact}) and the adjacency matrix \code{adj}.
#' @export
dag_pcts <- function(X, lags, alpha = 0.01, maxl = 3, maxc = 6) {
  X <- as.matrix(X); pc_ts_impl(X, lags, ncol(X), alpha = alpha, maxl = maxl, maxc = maxc)
}
