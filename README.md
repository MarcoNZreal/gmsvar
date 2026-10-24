# gmsvar

Graphical identification of structural vector autoregressions (SVARs) from the data: an extended-BIC search over contemporaneous orderings, the equivalence class under the time ordering, and the estimated SVAR. An optional first stage restricts the search to a sparse conditional independence graph (CIG) of the window of current and lagged variables, to impose prior knowledge or save computation.

```r
# install.packages("remotes"); remotes::install_github("MarcoNZreal/gmsvar")
library(gmsvar)
set.seed(1); X <- simulate_svar(svar_example(), 600); colnames(X) <- c("a", "b", "c")

fit <- gmsvar(X, pmax = 5)                        # lag order (AIC) -> extended-BIC search over orderings -> estimated SVAR
# optional: restrict the search with a CIG, your own (k(p+1) x k(p+1), or k x k) or an estimated one
E   <- cig_estimate(X, p = 2, method = "GMTS-ua", alpha = 0.3)
fit <- gmsvar(X, p = 2, cig = E)
summary(fit)        # final DAG: edges, coefficients, standard errors, compelled/reversible orientation, shock s.d.
svar_edges(fit); coef(fit); svar_irf(fit, 12)    # edge table, B0/A/sig2, impulse responses
plot(fit); svar_dot(fit)                          # picture, Graphviz source
gmsvar_boot(fit, X, B = 200)                      # block-bootstrap stability of every edge and orientation
```

* `var_order()`: AIC, BIC and Hannan-Quinn lag-order selection on a common sample. Checked against `vars::VARselect` (same selections in 90 of 90 test series; criteria differ by at most 1e-12 up to an additive constant).
* `cig_estimate()`: optional initial CIG of the window (test rule or penalised screen), to be passed as `cig =`.
* `cig_stats()`, `cig_select()`: window partial correlations and edge-selection rules (t-test, Bonferroni, SIN-type Sidak-Holm, BH, BIC cutoff). All rules rank the pairs identically and differ only in the cutoff.
* `cig_screen()`: adaptive graphical lasso and (adaptive) node-wise lasso screens.
* `dag_exact()`, `dag_heuristic()`: BIC or extended-BIC search (`ebic_gamma`) over contemporaneous orderings, optionally restricted to a CIG. `gmsvar()` uses the exact search up to k = 10 current variables and the heuristic above; the heuristic is not covered by the consistency result.
* `gmsvar()`: pipeline; returns the lag order, the initial CIG used, the DAG (`parents`), the estimated structural VAR (`svar`: `B0`, `A`, standard errors, shock variances) and the arrows that are common to the whole equivalence class (`compelled`). `svar_estimate()`, `svar_edges()`, `svar_irf()`, `svar_dot()`, `plot()`, `gmsvar_boot()` work on it.
* `dag_pcts()`: the tiered PC-stable comparator used in the paper (our own implementation, not tigramite).
* `random_svar()`, `simulate_svar()`, `true_parent_sets()`, `dag_compare()`: simulation and scoring.

Guidance from the accompanying paper: the default `gmsvar(X)` uses the extended BIC with `ebic_gamma = 0.5` (penalty per parent log N + 2 * 0.5 * log(m - 1)) and no screen. In the simulations (k = 3 to 50) it was about as accurate as the BIC for k <= 7 and far more accurate for k >= 20, where the BIC overfits (at k = 50, n = 1000: 13 against 72 skeleton errors), and it beat the adaptive node-wise lasso screen everywhere. Conservative test-based screens cost accuracy because the DAG step cannot restore dropped edges. Use `ebic_gamma = 0` for the plain BIC (version 0.3.0 behaviour with `screen = "auto"`). For very large systems a screen is a computational device: at k = 100 (m = 300) with N = 1000 the unscreened search was the most accurate but took about 50 minutes per fit, a liberal screen (`screen = "NW"`, or `cig_estimate(X, p, "GMTS-ua", alpha = 0.3)`) took 2-8 minutes with about twice the errors; with N close to m the node-wise lasso screen was about as accurate; with N < m only lasso screens can run and results are poor. A CIG is the way to impose prior knowledge.
