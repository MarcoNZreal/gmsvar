# dag.R -- automated CIG -> DAG step for SVARs, Markov-equivalence evaluation, Meek rules.
#
# Search: exact dynamic programming over contemporaneous orderings of the k current
# variables.  For each node i and each set S of allowed contemporaneous parents the
# parent set is chosen by backward elimination on BIC from the candidate set
# C = {CIG neighbours of i} within (S U lagged variables).  Lagged variables precede
# current ones (arrow of time).  Window indices: current 1..k, lagged k+1..m.

.rss <- function(Sg, i, P) {
  if (length(P) == 0) return(Sg[i, i])
  sol <- tryCatch(solve(Sg[P, P, drop = FALSE], Sg[P, i]), error = function(e) NULL)
  if (is.null(sol)) return(Inf)
  Sg[i, i] - sum(Sg[i, P] * sol)
}
.bic <- function(Sg, i, P, nv) {
  r <- .rss(Sg, i, P)
  if (!is.finite(r) || r <= 0) return(Inf)
  nv * log(r) + length(P) * log(nv)
}

best_parents <- function(Sg, i, cand, nv) {
  P <- cand
  cur <- .bic(Sg, i, P, nv)
  while (length(P) > 0) {
    best <- cur; bj <- NA
    for (j in P) {
      b <- .bic(Sg, i, setdiff(P, j), nv)
      if (b < best - 1e-12) { best <- b; bj <- j }
    }
    if (is.na(bj)) break
    P <- setdiff(P, bj); cur <- best
  }
  list(par = sort(P), score = cur)
}

# E: m x m logical candidate adjacency (CIG); E = NULL means no restriction.
# Returns list(parents = list over current nodes 1..k, score).
learn_dag <- function(X, lags, k, E = NULL, ebic_gamma = 0) {
  G <- window_mat(X, lags); nv <- nrow(G); m <- ncol(G)
  pen <- log(nv) + 2 * ebic_gamma * log(m - 1)   # BIC if ebic_gamma = 0, otherwise extended BIC with per-parent penalty
  Sg <- cov(G) * (nv - 1) / nv                    # ML covariance (ddof = 0)
  if (is.null(E)) E <- matrix(TRUE, m, m)
  E <- E | t(E)
  lagged <- if (m > k) (k + 1):m else integer(0)
  memo <- new.env(hash = TRUE)
  score <- function(i, allowed) {
    cand <- c(allowed, lagged); cand <- cand[cand != i & E[i, cand]]
    key <- paste0(i, ":", paste(cand, collapse = ","))
    v <- memo[[key]]
    if (is.null(v)) { v <- fast_best_parents(Sg, i, cand, nv, pen); memo[[key]] <- v }   # O(q^3) backward elimination by inverse downdating (best_parents() is the slow O(q^5) reference)
    v
  }
  full <- bitwShiftL(1L, k) - 1L
  bs <- numeric(full + 1); bi <- integer(full + 1); bp <- vector("list", full + 1)
  for (U in seq_len(full)) {
    members <- which(bitwAnd(U, bitwShiftL(1L, 0:(k - 1))) > 0)
    bsc <- Inf; bbi <- NA_integer_; bbp <- NULL
    for (i in members) {
      prev <- bitwAnd(U, bitwNot(bitwShiftL(1L, i - 1L)))
      allowed <- which(bitwAnd(prev, bitwShiftL(1L, 0:(k - 1))) > 0)
      sc <- score(i, allowed)
      tot <- bs[prev + 1] + sc$score
      if (tot < bsc - 1e-12) { bsc <- tot; bbi <- i; bbp <- sc$par }
    }
    bs[U + 1] <- bsc; bi[U + 1] <- bbi; bp[[U + 1]] <- bbp
  }
  parents <- vector("list", k); U <- full
  while (U > 0) {
    i <- bi[U + 1]; parents[[i]] <- bp[[U + 1]]
    U <- bitwAnd(U, bitwNot(bitwShiftL(1L, i - 1L)))
  }
  list(parents = parents, score = bs[full + 1])
}

true_parents <- function(m) {
  k <- nrow(m$A0)
  lapply(1:k, function(i) {
    P <- which(m$A0[i, ] != 0 & seq_len(k) != i)
    for (l in seq_along(m$A)) P <- c(P, l * k + which(m$A[[l]][i, ] != 0))
    sort(P)
  })
}

# ------------------------------------------------------------- evaluation
# Directed adjacency D[a,b] = TRUE if a -> b; lag-lag clique added (oriented by index)
.augmented <- function(par, k, m) {
  D <- matrix(FALSE, m, m)
  for (i in seq_len(k)) for (j in par[[i]]) D[j, i] <- TRUE
  if (m > k + 1) for (a in (k + 1):(m - 1)) for (b in (a + 1):m) D[a, b] <- TRUE
  D
}

.vstructs <- function(D) {
  adj <- D | t(D); out <- character(0)
  for (cc in seq_len(ncol(D))) {
    pa <- which(D[, cc])
    if (length(pa) >= 2) {
      cb <- combn(pa, 2)
      for (q in seq_len(ncol(cb))) {
        a <- cb[1, q]; b <- cb[2, q]
        if (!adj[a, b]) out <- c(out, paste(a, cc, b))
      }
    }
  }
  out
}

equivalent <- function(par1, par2, k, m) {
  D1 <- .augmented(par1, k, m); D2 <- .augmented(par2, k, m)
  isTRUE(all((D1 | t(D1)) == (D2 | t(D2)))) && setequal(.vstructs(D1), .vstructs(D2))
}

# Compelled / undirected status of contemporaneous edges via Meek rules R1-R3, with
# lag -> current edges known (time order).  Returns k x k logical 'dr' (a->b compelled).
cpdag_cc <- function(par, k, m) {
  D <- .augmented(par, k, m); adj <- D | t(D)
  dr <- D; und <- matrix(FALSE, m, m)
  for (a in seq_len(k)) for (b in seq_len(k)) if (a < b && adj[a, b]) {
    und[a, b] <- und[b, a] <- TRUE; dr[a, b] <- dr[b, a] <- FALSE
  }
  for (v in .vstructs(D)) {                         # a -> c <- b
    abc <- as.integer(strsplit(v, " ")[[1]]); a <- abc[1]; cc <- abc[2]; b <- abc[3]
    for (x in c(a, b)) if (und[x, cc]) { und[x, cc] <- und[cc, x] <- FALSE; dr[x, cc] <- TRUE }
  }
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (b in seq_len(k)) for (cc in seq_len(k)) {
      if (b == cc || !und[b, cc]) next
      # R1: a -> b, b - c, a not adjacent c
      r1 <- FALSE
      for (a in seq_len(m)) if (a != cc && dr[a, b] && !adj[a, cc]) { r1 <- TRUE; break }
      if (r1) { und[b, cc] <- und[cc, b] <- FALSE; dr[b, cc] <- TRUE; changed <- TRUE; next }
      # R2: b -> x -> c with b - c
      if (any(dr[b, ] & dr[, cc])) { und[b, cc] <- und[cc, b] <- FALSE; dr[b, cc] <- TRUE; changed <- TRUE; next }
      # R3: b - x1 -> c, b - x2 -> c, x1 not adjacent x2
      xs <- which(und[b, seq_len(k)] & dr[seq_len(k), cc])
      if (length(xs) >= 2) {
        cb <- combn(xs, 2)
        if (any(!adj[cbind(cb[1, ], cb[2, ])])) {
          und[b, cc] <- und[cc, b] <- FALSE; dr[b, cc] <- TRUE; changed <- TRUE
        }
      }
    }
  }
  dr[seq_len(k), seq_len(k), drop = FALSE]
}

dag_metrics <- function(par_est, par_true, k, m) {
  De <- .augmented(par_est, k, m); Dt <- .augmented(par_true, k, m)
  Ae <- De | t(De); At <- Dt | t(Dt)
  iu <- upper.tri(matrix(0, k, k))
  ec <- Ae[1:k, 1:k][iu]; tc <- At[1:k, 1:k][iu]
  el <- as.vector(Ae[1:k, (k + 1):m]); tl <- as.vector(At[1:k, (k + 1):m])
  out <- c(cc_TP = sum(ec & tc), cc_FP = sum(ec & !tc), cc_FN = sum(!ec & tc), cc_TN = sum(!ec & !tc),
           cl_TP = sum(el & tl), cl_FP = sum(el & !tl), cl_FN = sum(!el & tl), cl_TN = sum(!el & !tl))
  out["equiv"] <- as.integer(equivalent(par_est, par_true, k, m))
  out["skel_exact"] <- as.integer(out["cc_FP"] + out["cc_FN"] + out["cl_FP"] + out["cl_FN"] == 0)
  dt <- cpdag_cc(par_true, k, m); de <- cpdag_cc(par_est, k, m)
  idx <- which(dt, arr.ind = TRUE)
  right <- wrong <- undec <- miss <- 0
  for (q in seq_len(nrow(idx))) {
    a <- idx[q, 1]; b <- idx[q, 2]
    if (!Ae[a, b]) miss <- miss + 1
    else if (de[a, b]) right <- right + 1
    else if (de[b, a]) wrong <- wrong + 1
    else undec <- undec + 1
  }
  c(out, or_comp = nrow(idx), or_right = right, or_wrong = wrong, or_undec = undec, or_missing = miss)
}
