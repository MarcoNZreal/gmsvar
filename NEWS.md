# gmsvar 0.4.1

* `svar_estimate()` centres the window columns, so the reported coefficients come from the same regressions as the search (least squares with an intercept; standard errors with one degree of freedom for the mean). The change is in the fourth decimal for the application of the paper.
* The reference node routine `best_parents()`, a fallback of the fast routine, now uses the extended-BIC penalty when one is set.
* No change to the search: graphs, scores and timings are those of 0.4.0.

# gmsvar 0.4.0

* The default score of `gmsvar()` is now the extended BIC (`ebic_gamma = 0.5`, penalty per parent `log(N) + 2 * gamma * log(m - 1)`) and no screen (`screen = "none"`). In the accompanying simulations the BIC overfitted from about 15 variables, and the extended BIC was as accurate as, or more accurate than, the BIC with or without a screen at every size up to 50 variables. `ebic_gamma = 0` restores the BIC; `screen = "auto"` restores the screening rule of 0.3.0.
* `search = "auto"` now uses the full search over orderings (`dag_exact()`) up to 10 current variables (previously 7): with the faster node routine it takes about 1 s at k = 10 with two lags and 8 s with six, and it returned the same graph as the heuristic in 59 of 60 test fits at k = 8 to 10.
* `dag_exact()` and `dag_heuristic()` gain `ebic_gamma`; `gmsvar()` gains `restarts`; `dag_heuristic()` gains `restart = "perturb"` (iterated local search).
* The full search uses the O(q^3) backward-elimination node routine (identical parent sets to the previous routine, about 37 times faster); the heuristic search memoises node problems on the sorted predecessor set (about 4 times faster).
* When a variable has at least N - 1 candidate parents the node step switches to forward selection.
* Input checks: missing values, non-numeric data, constant or linearly dependent columns, invalid `ic`/`search`; an informative error when a partial-correlation screen is requested with more window variables than observations, and a warning when the unscreened search is run in that regime.
* Documentation: `dag_exact()` is described as a full search over orderings with backward-elimination parent selection, not an exact minimiser of the score.
* `gmsvar_boot()` carries the penalty and the number of restarts of the fit; `summary()` reports the criterion used.

# gmsvar 0.3.0

* Two-stage interface: `screen = "auto"` (no screen up to 20 variables, adaptive node-wise lasso above), `cig_estimate()`.

# gmsvar 0.2.0

* First public version: data-first `gmsvar()` pipeline, exact and heuristic BIC search over orderings, equivalence class under the time ordering, CIG tests and penalised screens, PC-ts, simulation tools, bootstrap, impulse responses.
