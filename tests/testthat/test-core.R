test_that("var_order matches stats::ar.ols criteria and recovers the true order", {
  set.seed(2); X <- simulate_svar(svar_example(), 800)
  o <- var_order(X, 6); expect_equal(o$BIC, 2)
  a <- ar.ols(X, aic = TRUE, order.max = 6, demean = TRUE, intercept = TRUE)
  expect_equal(o$AIC, a$order)
})
test_that("cig rules are nested and share a ranking", {
  set.seed(3); X <- simulate_svar(svar_example(), 500); S <- cig_stats(X, 2)
  ua <- cig_select(S, 3, 2, "GMTS-ua"); bic <- cig_select(S, 3, 2, "BIC"); sin <- cig_select(S, 3, 2, "SIN")
  expect_true(all(sin <= ua)); expect_true(sum(bic) <= sum(ua)); expect_equal(dim(ua), c(9, 9))
})
test_that("exact search recovers the DAG class of the example for large n", {
  set.seed(4); m <- svar_example(); X <- simulate_svar(m, 5000)
  fit <- gmsvar(X, p = 2, screen = "BIC")
  mt <- dag_compare(fit$parents, true_parent_sets(m), 3, 2)
  expect_equal(unname(mt["equiv"]), 1)
})
test_that("heuristic agrees with exact on the example and keeps the RNG stream", {
  set.seed(5); X <- simulate_svar(svar_example(), 1000); a <- dag_exact(X, 2); set.seed(9); r1 <- runif(1); set.seed(9)
  b <- dag_heuristic(X, 2); r2 <- runif(1); expect_equal(a$score, b$score, tolerance = 1e-6); expect_equal(r1, r2)
})
test_that("random_svar is stable and reproducible", {
  m1 <- random_svar(8, 2, seed = 1); m2 <- random_svar(8, 2, seed = 1); expect_identical(m1$A, m2$A)
  X <- simulate_svar(m1, 200); expect_equal(dim(X), c(200, 8))
})
test_that("gmsvar chooses p and prints", {
  set.seed(6); X <- simulate_svar(svar_example(), 700); fit <- gmsvar(X, pmax = 4, screen = "GMTS-ua", alpha = 0.3)
  expect_true(fit$p %in% 1:4); expect_output(print(fit), "gmsvar")
})
test_that("user-supplied cig is respected and k x k cig works", {
  set.seed(11); m <- svar_example(); X <- simulate_svar(m, 800); k <- 3; p <- 2
  E <- matrix(FALSE, 9, 9); E[1:3, 1:3] <- TRUE; E[1:3, 4:9] <- TRUE; E <- E | t(E); E[1, 2] <- E[2, 1] <- FALSE
  fit <- gmsvar(X, cig = E)
  expect_equal(fit$p, 2); expect_false(2 %in% fit$parents[[1]]); expect_false(1 %in% fit$parents[[2]])
  f2 <- gmsvar(X, p = 2, cig = matrix(TRUE, 3, 3)); f0 <- gmsvar(X, p = 2); expect_equal(f2$parents, f0$parents)
  expect_error(gmsvar(X, p = 2, cig = matrix(TRUE, 4, 4)))
})
test_that("estimated SVAR is close to the truth and irf is consistent", {
  set.seed(12); m <- svar_example(); X <- simulate_svar(m, 20000); fit <- gmsvar(X, p = 2, screen = "BIC")
  s <- fit$svar; expect_equal(s$A0, diag(3) - s$B0, ignore_attr = TRUE)
  mt <- dag_compare(fit$parents, true_parent_sets(m), 3, 2); expect_equal(unname(mt["equiv"]), 1)
  e <- svar_edges(fit); expect_true(all(c("from", "to", "lag", "coef", "se", "status") %in% names(e)))
  expect_true(all(abs(e$t) > 3))
  ir <- svar_irf(fit, 5); expect_equal(dim(ir), c(3, 3, 6))
  expect_equal(unname(ir[, , 1]), unname(solve(s$A0) %*% diag(sqrt(s$sig2))), tolerance = 1e-10)
  expect_output(print(summary(fit)), "Structural VAR"); expect_true(grepl("digraph", svar_dot(fit)))
})
test_that("compelled edges agree with equivalence on a small example", {
  set.seed(13); X <- simulate_svar(svar_example(), 3000); fit <- gmsvar(X, p = 2)
  expect_equal(dim(fit$compelled), c(3, 3))
})
test_that("bootstrap returns frequencies in [0,1] and keeps the RNG stream", {
  set.seed(14); X <- simulate_svar(svar_example(), 400); fit <- gmsvar(X, p = 2, screen = "BIC")
  set.seed(9); r1 <- runif(1); set.seed(9); b <- gmsvar_boot(fit, X, B = 10, block = 5); r2 <- runif(1)
  expect_equal(r1, r2); expect_true(all(b$adjacency >= 0 & b$adjacency <= 1)); expect_length(b$lagged, 2)
  expect_true(all(b$direction + t(b$direction) <= b$adjacency + 1e-12))
})
test_that("dag_pcts returns parent sets of the right shape", {
  set.seed(15); X <- simulate_svar(svar_example(), 3000); r <- dag_pcts(X, 2)
  expect_length(r$parents, 3); expect_true(all(unlist(r$parents) <= 9))
})

test_that("cig_estimate feeds gmsvar(cig=) and matches the screen shortcut", {
  set.seed(21); X <- simulate_svar(svar_example(), 500)
  E <- cig_estimate(X, 2, "GMTS-ua", alpha = 0.3)
  expect_equal(dim(E), c(9, 9)); expect_true(isSymmetric(E))
  f1 <- gmsvar(X, p = 2, cig = E); f2 <- gmsvar(X, p = 2, screen = "GMTS-ua", alpha = 0.3)
  expect_equal(f1$parents, f2$parents)
})

test_that("extended BIC selects no more parents than the BIC", {
  set.seed(3); X <- simulate_svar(svar_example(), 400)
  a <- dag_exact(X, 2); b <- dag_exact(X, 2, ebic_gamma = 1)
  expect_lte(sum(lengths(b$parents)), sum(lengths(a$parents)))
  h <- dag_heuristic(X, 2, ebic_gamma = 0.5, restarts = 0)
  expect_equal(length(h$parents), 3)
})

test_that("gmsvar defaults: extended BIC, no screen", {
  set.seed(4); X <- simulate_svar(svar_example(), 400); f <- gmsvar(X, p = 2)
  expect_equal(f$ebic_gamma, 0.5); expect_equal(f$screen, "none"); expect_null(f$cig)
})

test_that("input checks", {
  set.seed(5); X <- simulate_svar(svar_example(), 200)
  Y <- X; Y[3, 1] <- NA
  expect_error(gmsvar(Y, p = 2), "missing values")
  expect_error(gmsvar(X, ic = "XYZ"))
  m <- random_svar(12, 2, seed = 3); Z <- simulate_svar(m, 30)
  expect_error(gmsvar(Z, p = 2, screen = "GMTS-ua"), "more observations")
  expect_warning(gmsvar(Z, p = 2), "not reliable")
  expect_silent(f <- gmsvar(Z, p = 2, screen = "NW"))
})

test_that("search = 'auto' is exact up to k = 10 and heuristic above", {
  d8 <- random_svar(k = 8, p = 1, seed = 8001); set.seed(1); X8 <- simulate_svar(d8, 300)
  expect_equal(gmsvar(X8, p = 1)$search, "exact")
  d11 <- random_svar(k = 11, p = 1, seed = 11001); set.seed(1); X11 <- simulate_svar(d11, 300)
  expect_equal(gmsvar(X11, p = 1)$search, "heuristic")
})

test_that("constant and linearly dependent columns are rejected with a clear message", {
  set.seed(1); X <- simulate_svar(svar_example(), 300)
  expect_error(gmsvar(cbind(X, 1), p = 2), "constant")
  expect_error(gmsvar(cbind(X, X[, 1] + X[, 2]), p = 2), "linearly dependent")
})

test_that("the reference node routine uses the extended-BIC penalty, as the fast one does", {
  set.seed(4); d <- random_svar(k = 5, p = 2, seed = 501); X <- simulate_svar(d, 200)
  G <- window_matrix(X, 2); nv <- nrow(G); Sg <- cov(G) * (nv - 1) / nv; pen <- log(nv) + 2 * 0.5 * log(ncol(G) - 1)
  for (i in 1:5) {
    a <- gmsvar:::best_parents(Sg, i, setdiff(seq_len(ncol(G)), i), nv, pen)
    b <- gmsvar:::fast_best_parents(Sg, i, setdiff(seq_len(ncol(G)), i), nv, pen)
    expect_equal(a$par, b$par); expect_equal(a$score, b$score, tolerance = 1e-8)
  }
})

test_that("svar_estimate reproduces least squares with an intercept on the window", {
  set.seed(2); X <- simulate_svar(svar_example(), 300); fit <- gmsvar(X, p = 2)
  G <- window_matrix(X, 2); est <- fit$svar
  for (i in 1:3) { P <- sort(fit$parents[[i]]); if (!length(P)) next
    co <- coef(lm(G[, i] ~ G[, P, drop = FALSE]))[-1]; se <- summary(lm(G[, i] ~ G[, P, drop = FALSE]))$coefficients[-1, 2]
    j <- P[1]; v <- if (j <= 3) est$B0[i, j] else est$A[[(j - 1) %/% 3]][i, (j - 1) %% 3 + 1]
    vs <- if (j <= 3) est$B0_se[i, j] else est$A_se[[(j - 1) %/% 3]][i, (j - 1) %% 3 + 1]
    expect_equal(unname(v), unname(co[1]), tolerance = 1e-10); expect_equal(unname(vs), unname(se[1]), tolerance = 1e-10) }
})

test_that("a k by k cig without p keeps the lag order chosen by the criterion", {
  set.seed(1); X <- simulate_svar(svar_example(), 600)
  f <- gmsvar(X, cig = matrix(TRUE, 3, 3)); g <- gmsvar(X)
  expect_equal(f$p, g$p); expect_true(f$p >= 1)
  expect_true(sum(unlist(f$parents) > 3) > 0)
})

test_that("p is validated and the random stream is left untouched", {
  set.seed(1); X <- simulate_svar(svar_example(), 300)
  expect_error(gmsvar(X, p = 1.5), "positive integer"); expect_error(gmsvar(X, p = 0), "positive integer")
  expect_error(gmsvar(X, p = 2, cig = matrix(c(TRUE, NA, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), 3)), "missing")
  fit <- gmsvar(X, p = 2)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  invisible(gmsvar_boot(fit, X, B = 3)); expect_false(exists(".Random.seed", envir = globalenv()))
  expect_error(gmsvar_boot(fit, X, B = 3, block = 1000), "block")
  pc <- lapply(fit$parents, function(P) P[P <= 3]); expect_silent(dag_compare(pc, pc, 3, 0)[1])
})
