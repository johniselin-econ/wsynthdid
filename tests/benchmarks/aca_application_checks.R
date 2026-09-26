# =============================================================================
# aca_application_checks.R — weighted SDID checks on the ACA application data
# (formerly tests/test_weighted.R). Not part of R CMD check: it needs
# paper/data/analysis_data.csv from the paper repository
# (johniselin-econ/weighted-sdid), which left this repo at the July 2026 split.
# To run: copy that file to ../paper/data/ relative to this directory, then
#   Rscript aca_application_checks.R
# Last recorded output: aca_application_checks_output.txt
# =============================================================================

# ── Load package functions ────────────────────────────────────────────────────
invisible(lapply(list.files("../R", pattern = "[.][Rr]$", full.names = TRUE), source))

# ── Load and prepare data (mirrors analyze_application.R) ────────────────────
panel <- read.csv("../paper/data/analysis_data.csv")
panel <- panel[, c("fips", "year", "crude_rate", "expansion", "population")]

# Build panel matrices (same logic as analyze_application.R)
dfp <- data.frame(
  .unit = as.factor(panel$fips),
  .time = panel$year,
  y     = panel$crude_rate,
  .W    = as.integer(panel$expansion == 1 & panel$year >= 2014)
)
setup <- panel.matrices(dfp)
Y     <- setup$Y
N0    <- setup$N0
T0    <- setup$T0
N1    <- nrow(Y) - N0

# Population weights (2013, same as paper)
pop_2013      <- panel[panel$expansion == 1 & panel$year == 2013, c("fips", "population")]
treated_names <- rownames(Y)[(N0 + 1):nrow(Y)]
pop_map       <- setNames(pop_2013$population, as.character(pop_2013$fips))
treated_weights <- as.numeric(pop_map[treated_names])
treated_weights <- treated_weights / sum(treated_weights, na.rm = TRUE)
uniform_weights <- rep(1 / N1, N1)

cat("Data loaded:", nrow(Y), "units,", ncol(Y), "periods,", N1, "treated\n\n")

# =============================================================================
# 1. Unweighted baseline
# =============================================================================
tau_sdid = synthdid_estimate(Y, N0, T0)
stopifnot(is.finite(as.numeric(tau_sdid)))
cat("PASS 1: unweighted ATT finite:", round(as.numeric(tau_sdid), 3), "\n")

# =============================================================================
# 2. Weighted estimate is finite
# =============================================================================
tau_w = synthdid_estimate_weighted(Y, N0, T0, treated.weights = treated_weights)
stopifnot(is.finite(as.numeric(tau_w)))
cat("PASS 2: weighted ATT finite:", round(as.numeric(tau_w), 3), "\n")

# =============================================================================
# 3. Uniform treated weights recover unweighted estimate
# =============================================================================
tau_w_unif = synthdid_estimate_weighted(Y, N0, T0, treated.weights = uniform_weights)
stopifnot(abs(as.numeric(tau_w_unif) - as.numeric(tau_sdid)) < 1e-4)
cat("PASS 3: uniform weighted == unweighted\n")

# =============================================================================
# 4. Population weights produce a different ATT than uniform weights
#    (paper result: tau_sdid ~ 0, tau_w ~ -17.5)
# =============================================================================
stopifnot(abs(as.numeric(tau_w) - as.numeric(tau_sdid)) > 1.0)
cat("PASS 4: population weighted differs from unweighted by >1\n")

# =============================================================================
# 5. Omega and lambda sum to one
# =============================================================================
w = attr(tau_w, 'weights')
stopifnot(abs(sum(w$omega)  - 1) < 1e-6)
stopifnot(abs(sum(w$lambda) - 1) < 1e-6)
cat("PASS 5: omega and lambda sum to one\n")

# =============================================================================
# 6. format() produces a non-empty string (tests our fix #1)
# =============================================================================
out = format(tau_w)
stopifnot(is.character(out), nchar(out) > 0)
cat("PASS 6: format.synthdid_estimate_weighted works\n")

# =============================================================================
# 7. summary() returns expected fields
# =============================================================================
s = summary(tau_w, fast = TRUE)
stopifnot(all(c('estimate', 'se', 'controls', 'periods', 'dimensions') %in% names(s)))
stopifnot(is.finite(s$estimate))
cat("PASS 7: summary fields present\n")

# =============================================================================
# 8. Unweighted jackknife SE is positive and finite
# =============================================================================
se_jk_uw = sqrt(vcov(tau_sdid, method = 'jackknife'))
stopifnot(is.finite(se_jk_uw), se_jk_uw > 0)
cat("PASS 8: unweighted jackknife SE positive:", round(se_jk_uw, 3), "\n")

# =============================================================================
# 8a. plot() dispatch works (tests fix #4)
# =============================================================================
p = plot(tau_w)
stopifnot(!is.null(p))
cat("PASS 8a: plot.synthdid_estimate_weighted works\n")

# =============================================================================
# 8b. Event study returns correct shape (tests fix #4 S3 dispatch)
# =============================================================================
es = synthdid_event_study(tau_w)
stopifnot(inherits(es, 'synthdid_event_study'))
stopifnot(all(c('relative_time', 'estimate') %in% names(es)))
tryCatch(plot(es), error = function(e) stop("plot.synthdid_event_study failed: ", e$message))
cat("PASS 8b: event study shape and plot dispatch work\n")

# =============================================================================
# 9. Jackknife SE is positive and finite (requires N1 > 1 — passes here)
# =============================================================================
se_jk = sqrt(vcov(tau_w, method = 'jackknife'))
stopifnot(is.finite(se_jk), se_jk > 0)
cat("PASS 9: jackknife SE positive:", round(se_jk, 3), "\n")

# =============================================================================
# 10. Bootstrap SE is positive and finite
# =============================================================================
se_bt = sqrt(vcov(tau_w, method = 'bootstrap'))
stopifnot(is.finite(se_bt), se_bt > 0)
cat("PASS 10: bootstrap SE positive:", round(se_bt, 3), "\n")

# =============================================================================
# 11. Placebo SE is positive and finite
# =============================================================================
se_pl = sqrt(vcov(tau_w, method = 'placebo'))
stopifnot(is.finite(se_pl), se_pl > 0)
cat("PASS 11: placebo SE positive:", round(se_pl, 3), "\n")

# =============================================================================
# 12. DID weighted estimate is finite
# =============================================================================
tau_did_w = did_estimate_weighted(Y, N0, T0, treated.weights = treated_weights)
stopifnot(is.finite(as.numeric(tau_did_w)))
cat("PASS 12: did_estimate_weighted finite:", round(as.numeric(tau_did_w), 3), "\n")

# =============================================================================
# 13. SC weighted estimate is finite
# =============================================================================
tau_sc_w = sc_estimate_weighted(Y, N0, T0, treated.weights = treated_weights)
stopifnot(is.finite(as.numeric(tau_sc_w)))
cat("PASS 13: sc_estimate_weighted finite:", round(as.numeric(tau_sc_w), 3), "\n")

# =============================================================================
# 15. Effective N1 is less than raw N1 (sanity check on weight concentration)
# =============================================================================
n1_eff = 1 / sum(treated_weights^2)
stopifnot(n1_eff < N1)
cat("PASS 14: N1_eff", round(n1_eff, 0), "< N1", N1, "\n")

stopifnot(as.numeric(tau_w) > -25 & as.numeric(tau_w) < -5)
cat("PASS 15: weighted ATT", round(as.numeric(tau_w), 3), "in expected range (-25, -5)\n")

cat("\n=== All 15 tests passed ===\n")
