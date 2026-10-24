# dag_heur.R -- heuristic CIG -> DAG search for larger k: local search over contemporaneous orderings.
# Node problems are the same as in the exact DP (backward elimination on BIC from the candidate set, memoised on (i, candidate set));
# backward elimination uses inverse downdating (O(q^2) per removal) and gives the same answer as best_parents() in dag.R.

fast_best_parents <- function(Sg, i, cand, nv, pen = log(nv)) {
  P <- cand; q <- length(P); lg <- pen   # penalty per parent: log N (BIC) or log N + 2 gamma log(m - 1) (extended BIC)
  if (q == 0) return(list(par = integer(0), score = nv * log(Sg[i, i])))
  if (q >= nv - 1) return(forward_parents(Sg, i, cand, nv, pen))   # more candidates than observations: the full regression is not estimable, select forwards
  Om <- tryCatch(chol2inv(chol(Sg[P, P, drop = FALSE])), error = function(e) NULL)
  if (is.null(Om)) return(best_parents(Sg, i, cand, nv))
  s <- Sg[P, i]; rss <- Sg[i, i] - sum(s * (Om %*% s))
  if (!is.finite(rss) || rss <= 0) return(list(par = sort(P), score = Inf))
  cur <- nv * log(rss) + q * lg
  repeat {
    if (length(P) == 0) break
    beta <- as.vector(Om %*% s); d <- diag(Om)
    inc <- beta^2 / d                         # increase of RSS when dropping each variable
    jb <- which.min(inc); rss_new <- rss + inc[jb]
    nb <- nv * log(rss_new) + (length(P) - 1) * lg
    if (!is.finite(nb)) break
    if (nb < cur - 1e-12) {
      oj <- Om[, jb]; Om <- Om[-jb, -jb, drop = FALSE] - tcrossprod(oj[-jb]) / oj[jb]
      P <- P[-jb]; s <- s[-jb]; rss <- rss_new; cur <- nb
    } else break
  }
  list(par = sort(P), score = cur)
}

learn_dag_heur <- function(X, lags, k, E = NULL, restarts = 3, seed = 1, maxpass = 50, restart = c("random", "perturb"), strength = max(3L, round(k / 6)), ebic_gamma = 0) {
  restart <- match.arg(restart)
  G <- window_mat(X, lags); nv <- nrow(G); m <- ncol(G)
  pen <- log(nv) + 2 * ebic_gamma * log(m - 1)   # extended BIC (Chen and Chen 2008) with the binomial term approximated per parent, so the penalty does not depend on the ordering
  Sg <- cov(G) * (nv - 1) / nv
  if (is.null(E)) E <- matrix(TRUE, m, m)
  E <- E | t(E)
  lagged <- if (m > k) (k + 1):m else integer(0)
  memo <- new.env(hash = TRUE); nev <- 0L; nhit <- 0L
  node <- function(i, allowed) {
    cand <- c(sort(allowed), lagged); cand <- cand[cand != i & E[i, cand]]   # sorted: the same set reached through different orders hits the same memo entry
    key <- paste0(i, ":", paste(cand, collapse = ","))
    v <- memo[[key]]
    if (is.null(v)) { nev <<- nev + 1L; v <- fast_best_parents(Sg, i, cand, nv, pen); memo[[key]] <- v } else nhit <<- nhit + 1L
    v
  }
  ord_score <- function(ord) { sc <- numeric(k); for (r in seq_len(k)) sc[ord[r]] <- node(ord[r], ord[seq_len(r - 1)])$score; sc }
  greedy <- function() {
    placed <- integer(0); rest <- seq_len(k)
    base <- vapply(seq_len(k), function(i) node(i, integer(0))$score, 0)
    while (length(rest)) {
      gain <- vapply(rest, function(i) base[i] - node(i, placed)$score, 0)
      j <- rest[which.min(gain)]; placed <- c(placed, j); rest <- setdiff(rest, j)
    }
    placed
  }
  improve <- function(ord) {
    # moves: insertion of one node at any other position (includes adjacent swaps); first-improvement until local optimum
    sc <- ord_score(ord); cur <- sum(sc)
    for (pass in seq_len(maxpass)) {
      improved <- FALSE
      for (a in sample(ord)) {
        r0 <- match(a, ord); rest <- ord[-r0]
        bestv <- cur; bestord <- NULL; bestsc <- NULL
        for (r in seq_len(k)) {
          if (r == r0) next
          o2 <- append(rest, a, after = r - 1)
          lo <- min(r, r0); hi <- max(r, r0); s2 <- sc
          for (pos in lo:hi) s2[o2[pos]] <- node(o2[pos], o2[seq_len(pos - 1)])$score
          if (sum(s2) < bestv - 1e-9) { bestv <- sum(s2); bestord <- o2; bestsc <- s2 }
        }
        if (!is.null(bestord)) { ord <- bestord; sc <- bestsc; cur <- bestv; improved <- TRUE }
      }
      if (!improved) break
    }
    list(ord = ord, score = cur)
  }
  # use a private RNG stream: save the caller's .Random.seed and restore it on exit (set.seed here used to reset the caller's stream)
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  best <- improve(greedy())
  for (x in seq_len(restarts)) {
    if (restart == "random") o <- sample.int(k)
    else { o <- best$ord; for (j in seq_len(strength)) { a <- sample.int(k, 1); o <- append(o[-match(a, o)], a, after = sample.int(k, 1) - 1) } }   # iterated local search: perturb the incumbent by random insertions
    r <- improve(o); if (r$score < best$score) best <- r
  }
  parents <- vector("list", k)
  for (r in seq_len(k)) parents[[best$ord[r]]] <- node(best$ord[r], best$ord[seq_len(r - 1)])$par
  list(parents = parents, score = best$score, order = best$ord, evals = nev, hits = nhit)
}

# Forward selection on the (extended) BIC, used when the candidate set is too large for backward elimination (q >= N - 1).
forward_parents <- function(Sg, i, cand, nv, pen = log(nv)) {
  P <- integer(0); cur <- nv * log(Sg[i, i]); rest <- cand
  repeat {
    if (length(rest) == 0 || length(P) >= nv - 2) break
    sc <- vapply(rest, function(j) { Q <- c(P, j); r <- tryCatch(Sg[i, i] - sum(Sg[Q, i] * solve(Sg[Q, Q, drop = FALSE], Sg[Q, i])), error = function(e) NA_real_)
      if (!is.finite(r) || r <= 0) Inf else nv * log(r) + length(Q) * pen }, 0)
    jb <- which.min(sc); if (sc[jb] >= cur - 1e-12) break
    P <- c(P, rest[jb]); rest <- rest[-jb]; cur <- sc[jb]
  }
  list(par = sort(P), score = cur)
}
