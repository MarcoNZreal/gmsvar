# gmsvar

Data-driven identification of structural vector autoregressions (SVARs): lag order, BIC search over contemporaneous DAGs, equivalence class under the time ordering, estimated SVAR. A conditional independence graph (CIG) is an optional input.

```r
# install.packages("remotes"); remotes::install_github("<user>/gmsvar")   # once published
library(gmsvar)
set.seed(1); X <- simulate_svar(svar_example(), 600); colnames(X) <- c("a", "b", "c")

fit <- gmsvar(X, pmax = 5)                        # lag order (AIC) -> BIC search over orderings -> estimated SVAR
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
* `dag_exact()`, `dag_heuristic()`: BIC search over contemporaneous orderings, optionally restricted to a CIG. The heuristic is not covered by the consistency result.
* `gmsvar()`: pipeline; returns the lag order, the initial CIG used, the DAG (`parents`), the estimated structural VAR (`svar`: `B0`, `A`, standard errors, shock variances) and the arrows that are common to the whole equivalence class (`compelled`). `svar_estimate()`, `svar_edges()`, `svar_irf()`, `svar_dot()`, `plot()`, `gmsvar_boot()` work on it.
* `dag_pcts()`: the tiered PC-stable comparator used in the paper (our own implementation, not tigramite).
* `random_svar()`, `simulate_svar()`, `true_parent_sets()`, `dag_compare()`: simulation and scoring.

Guidance from the accompanying paper: the unrestricted search from the data was the most accurate; a CIG helps mainly to impose prior knowledge or, with a liberal screen, to speed up the search for large k. Conservative screens cost accuracy because the DAG step cannot restore dropped edges.
