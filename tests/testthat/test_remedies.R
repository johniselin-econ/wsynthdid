# =============================================================================
# Tests for the trend-robust variants: detrended and stratified weighted SDID
# =============================================================================

# A small panel with unit-specific linear trends whose slopes are correlated
# with the treated-unit weights, plus a constant treatment effect tau. With
# treated.slope = 0 the treated trend is inside the convex hull of control
# slopes and SDID matches it, so estimates are unbiased and the detrended /
# stratified variants should agree with the raw fit. treated.slope adds a
# common differential trend to the treated units; it biases comparison-based
# estimates (by roughly the unmatched slope x post-minus-pre time gap) only
# when it pushes the weighted treated slope OUTSIDE the control-slope hull —
# with trend.scale = 0.5 the size-correlated control slopes have wide support,
# so use trend.scale = 0 (all control slopes ~0) to make any positive
# treated.slope genuinely unmatchable.
trend.dgp = function(N0 = 40, N1 = 12, T0 = 8, T1 = 4, tau = 1, trend.scale = 0.5,
                     treated.slope = 0, sigma = 0.02, seed = 42) {
  set.seed(seed)
  N = N0 + N1; T = T0 + T1
  sizes  = exp(rnorm(N, sd = 1))                      # unit sizes, all units
  slopes = trend.scale * (sizes / mean(sizes) - 1) +  # size-correlated trends
    treated.slope * (seq_len(N) > N0)                 # + differential treated trend
  level  = rnorm(N)
  W = matrix(0, N, T); W[(N0 + 1):N, (T0 + 1):T] = 1
  Y = outer(level, rep(1, T)) + outer(slopes, 1:T) + tau * W +
    matrix(rnorm(N * T, sd = sigma), N, T)
  rownames(Y) = 1:N; colnames(Y) = 1:T
  tw = sizes[(N0 + 1):N] / sum(sizes[(N0 + 1):N])
  list(Y = Y, N0 = N0, T0 = T0, tw = tw, sizes = sizes, tau = tau)
}

test_that("unit_linear_trend reproduces an exactly linear panel", {
  N = 10; T = 9; T0 = 5
  Y = outer(rnorm(N), rep(1, T)) + outer(rnorm(N), 1:T)
  expect_equal(unit_linear_trend(Y, T0), Y, tolerance = 1e-10)
})

test_that("detrend = TRUE equals estimating on manually detrended data", {
  d = trend.dgp()
  tau.detrend = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw,
                                           detrend = TRUE)
  Yd = d$Y - unit_linear_trend(d$Y, d$T0)
  tau.manual = synthdid_estimate_weighted(Yd, d$N0, d$T0, treated.weights = d$tw)
  expect_equal(c(tau.detrend), c(tau.manual), tolerance = 1e-10)
  # setup keeps the RAW outcome; detrend is recorded in the refit options
  expect_equal(attr(tau.detrend, 'setup')$Y, d$Y)
  expect_true(isTRUE(attr(tau.detrend, 'opts')$detrend))
})

test_that("detrending removes unmatchable differential-trend bias", {
  d = trend.dgp(trend.scale = 0, treated.slope = 0.3)
  tau.raw = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw)
  tau.det = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw,
                                       detrend = TRUE)
  expect_gt(abs(c(tau.raw) - d$tau), 0.5)          # raw estimate is visibly biased
  expect_lt(abs(c(tau.det) - d$tau), 0.05)         # detrended recovers tau
})

test_that("detrend is (approximately) harmless on trend-free data", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = rep(1 / N1, N1)
  tau.raw = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  tau.det = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw,
                                       detrend = TRUE)
  expect_lt(abs(c(tau.raw) - c(tau.det)), 0.25 * abs(c(tau.raw)) + 0.1)
})

test_that("detrend errors with covariates or too-short pre-period", {
  d = trend.dgp()
  X = array(rnorm(length(d$Y)), dim = c(dim(d$Y), 1))
  expect_error(synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw,
                                          X = X, detrend = TRUE),
               "covariates")
  expect_error(synthdid_estimate_weighted(d$Y[, 1:6], d$N0, 2, treated.weights = d$tw,
                                          detrend = TRUE),
               "T0 >= 3")
})

test_that("effect curve of a detrended estimate averages back to the estimate", {
  d = trend.dgp()
  tau.det = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw,
                                       detrend = TRUE)
  curve = synthdid_effect_curve_weighted(tau.det)
  pw = attr(tau.det, 'period.weights')
  expect_equal(sum(curve * pw), c(tau.det), tolerance = 1e-10)
})

test_that("bootstrap vcov of a detrended estimate re-fits slopes and runs", {
  d = trend.dgp()
  cl = rep(1:13, length.out = nrow(d$Y))
  tau.det = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw,
                                       cluster = cl, detrend = TRUE)
  set.seed(1)
  v = vcov(tau.det, method = "bootstrap", replications = 20)
  expect_true(is.finite(v[1]) && v[1] > 0)
})

test_that("stratified estimate with a single stratum equals the weighted estimate", {
  d = trend.dgp()
  tau.w = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw)
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = rep(1, nrow(d$Y)),
                                       treated.weights = d$tw)
  expect_equal(c(tau.w), as.numeric(tau.s), tolerance = 1e-10)
})

test_that("stratified aggregate is the share-weighted sum of stratum estimates", {
  d = trend.dgp()
  strata = cut(d$sizes, quantile(d$sizes, c(0, .5, 1)), include.lowest = TRUE, labels = FALSE)
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw)
  tab = attr(tau.s, 'strata.table')
  expect_equal(sum(tab$share), 1, tolerance = 1e-10)
  expect_equal(as.numeric(tau.s), sum(tab$share * tab$estimate), tolerance = 1e-10)
  # shares are the treated-weight mass of each stratum (same estimand)
  strata.treated = strata[(d$N0 + 1):nrow(d$Y)]
  for (s in tab$stratum) {
    expect_equal(tab$share[tab$stratum == s],
                 sum(d$tw[as.character(strata.treated) == s]), tolerance = 1e-10)
  }
})

test_that("stratified estimate recovers tau under size-correlated trends", {
  d = trend.dgp()
  strata = cut(d$sizes, quantile(d$sizes, seq(0, 1, .25)), include.lowest = TRUE, labels = FALSE)
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw)
  expect_lt(abs(as.numeric(tau.s) - d$tau), 0.1)
})

test_that("infeasible stratum errors, or is dropped with drop.infeasible", {
  d = trend.dgp()
  strata = rep(1, nrow(d$Y))
  strata[nrow(d$Y)] = 2                 # last treated unit alone: no controls in stratum 2
  expect_error(synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                            treated.weights = d$tw),
               "stratum")
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw, drop.infeasible = TRUE)
  expect_equal(attr(tau.s, 'strata.dropped'), "2")
  expect_equal(sum(attr(tau.s, 'strata.table')$share), 1, tolerance = 1e-10)
})

test_that("cluster bootstrap vcov runs for the stratified estimator", {
  d = trend.dgp()
  strata = cut(d$sizes, quantile(d$sizes, c(0, .5, 1)), include.lowest = TRUE, labels = FALSE)
  cl = rep(1:13, length.out = nrow(d$Y))
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw, cluster = cl)
  set.seed(1)
  v = vcov(tau.s, replications = 20)
  expect_true(is.finite(v[1]) && v[1] > 0)
})

test_that("stratified passes detrend through to within-stratum fits", {
  d = trend.dgp()
  strata = cut(d$sizes, quantile(d$sizes, c(0, .5, 1)), include.lowest = TRUE, labels = FALSE)
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw, detrend = TRUE)
  fits = attr(tau.s, 'strata.fits')
  expect_true(all(vapply(fits, function(f) isTRUE(attr(f, 'opts')$detrend), logical(1))))
  expect_lt(abs(as.numeric(tau.s) - d$tau), 0.1)
})

test_that("synthdid_event_study can return the replication draws", {
  d = trend.dgp()
  tau.w = synthdid_estimate_weighted(d$Y, d$N0, d$T0, treated.weights = d$tw)
  set.seed(1)
  es = synthdid_event_study(tau.w, se.method = "bootstrap", replications = 10,
                            return.replications = TRUE)
  draws = attr(es, 'replications')
  expect_equal(dim(draws), c(10, ncol(d$Y)))
  expect_equal(apply(draws, 2, sd), es$se, tolerance = 1e-10)
})

test_that("stratified fit reports per-stratum n1_eff and dropped treated mass", {
  d = trend.dgp()
  strata = cut(d$sizes, quantile(d$sizes, c(0, .5, 1)), include.lowest = TRUE, labels = FALSE)
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw)
  tab = attr(tau.s, 'strata.table')
  strata.treated = strata[(d$N0 + 1):nrow(d$Y)]
  for (s in tab$stratum) {
    w = d$tw[as.character(strata.treated) == s]; w = w / sum(w)
    expect_equal(tab$n1_eff[tab$stratum == s], 1 / sum(w^2), tolerance = 1e-10)
  }
  expect_equal(attr(tau.s, 'dropped.mass'), 0, tolerance = 1e-12)

  strata2 = rep(1, nrow(d$Y)); strata2[nrow(d$Y)] = 2
  tau.d = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata2,
                                       treated.weights = d$tw, drop.infeasible = TRUE)
  tw.norm = d$tw / sum(d$tw)
  expect_equal(attr(tau.d, 'dropped.mass'), tw.norm[length(tw.norm)], tolerance = 1e-10)
})

test_that("stratified_boot_rep return modes agree", {
  d = trend.dgp()
  strata = cut(d$sizes, quantile(d$sizes, c(0, .5, 1)), include.lowest = TRUE, labels = FALSE)
  cl = rep(1:10, length.out = nrow(d$Y))
  tau.s = synthdid_estimate_stratified(d$Y, d$N0, d$T0, strata = strata,
                                       treated.weights = d$tw, cluster = cl)
  set.seed(7); e = stratified_boot_rep(tau.s, "cluster", cl)
  set.seed(7); dt = stratified_boot_rep(tau.s, "cluster", cl, return = "detail")
  set.seed(7); ob = stratified_boot_rep(tau.s, "cluster", cl, return = "object")
  expect_equal(unname(dt["estimate"]), e)
  expect_true(dt["dropped.mass"] >= 0 && dt["dropped.mass"] < 1)
  expect_s3_class(ob, "synthdid_estimate_stratified")
  expect_equal(as.numeric(ob), e)
})

test_that("weighted omega minimizes the weighted unit-level pre-fit criterion (Prop 0')", {
  # With the penalty fixed, the collapsed-target objective and the tw-weighted
  # unit-level objective differ by a constant in omega, so the two programs have
  # the same solution. Check via the objective difference at two feasible points.
  d = trend.dgp()
  N0 = d$N0; T0 = d$T0; N1 = nrow(d$Y) - N0
  tw = d$tw / sum(d$tw)
  Yco = d$Y[1:N0, 1:T0]; Ytr = d$Y[(N0 + 1):nrow(d$Y), 1:T0, drop = FALSE]
  target = colSums(tw * Ytr)
  f_coll = function(w0, w) sum((w0 + colSums(w * Yco) - target)^2)
  f_unit = function(w0, w) sum(tw * rowSums(sweep(-Ytr, 2, w0 + colSums(w * Yco), "+")^2))
  set.seed(3)
  a = runif(N0); a = a / sum(a); b = runif(N0); b = b / sum(b)
  expect_equal(f_unit(0.3, a) - f_coll(0.3, a), f_unit(-1, b) - f_coll(-1, b), tolerance = 1e-8)
})
