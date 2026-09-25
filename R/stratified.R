# =============================================================================
# STRATIFIED WEIGHTED SDID
# Weighted SDID within researcher-defined strata (donor pool restricted to
# same-stratum controls), aggregated with treated-weight stratum shares.
# The aggregate targets the SAME estimand tau(omega-tilde) as the unrestricted
# weighted estimator: only the comparisons are restricted, in the spirit of
# Abadie & L'Hour (2021) penalized SC and Abadie (2021) donor-pool advice.
# =============================================================================

#' Stratified weighted SDID estimate.
#'
#' Partitions units into strata (e.g. county-size bins), runs weighted SDID
#' within each stratum with the donor pool restricted to same-stratum controls,
#' and aggregates the within-stratum estimates using the treated-weight share
#' of each stratum. Because the aggregation weights are the researcher-chosen
#' treated-unit weights, the aggregate estimates the same weighted estimand as
#' [synthdid_estimate_weighted]; stratification only restricts which controls
#' each treated unit may be compared to.
#'
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param strata a vector of length nrow(Y) (controls first, like Y) giving each unit's stratum.
#'        Coerced to factor. Every stratum that contains a treated unit must also contain at
#'        least `min.controls` control units, unless drop.infeasible = TRUE.
#' @param treated.weights A vector of weights for treated units (length N1 = nrow(Y) - N0).
#'        If NULL, uses uniform weights. Normalized to sum to 1.
#' @param cluster An optional vector of length nrow(Y) of cluster IDs (e.g. state FIPS).
#'        Stored for cluster-robust vcov(); also subset into each within-stratum fit.
#' @param min.controls Minimum number of same-stratum controls required for a stratum containing
#'        treated units to be estimable. Default 2.
#' @param drop.infeasible If TRUE, strata with treated units but fewer than `min.controls`
#'        controls (or no treated units after resampling) are dropped and the remaining stratum
#'        shares renormalized, rather than raising an error. Used by the bootstrap, where a
#'        resample can empty one side of a stratum. Default FALSE.
#' @param ... additional options passed to every within-stratum [synthdid_estimate_weighted]
#'        call (e.g. detrend = TRUE).
#' @return The aggregate estimate: a scalar of class 'synthdid_estimate_stratified' with attributes
#'         'strata.fits' (named list of within-stratum synthdid_estimate_weighted objects),
#'         'strata.table' (data.frame: stratum, N1, N0, share, estimate, n1_eff --
#'         the Kish effective number of treated units within the stratum),
#'         'strata.dropped' (names of dropped strata) and 'dropped.mass' (the
#'         treated-weight mass they carried; > 0 means the aggregate targets the
#'         renormalized rather than the original estimand),
#'         'setup' (list: Y, N0, T0, strata), 'treated.weights', 'cluster', and 'opts'
#'         (the extra arguments, re-applied on every vcov() refit).
#' @export synthdid_estimate_stratified
synthdid_estimate_stratified = function(Y, N0, T0, strata,
                                        treated.weights = NULL,
                                        cluster = NULL,
                                        min.controls = 2,
                                        drop.infeasible = FALSE,
                                        ...) {
  N  = nrow(Y)
  N1 = N - N0
  stopifnot(N1 >= 1, length(strata) == N, is.null(cluster) || length(cluster) == N)
  strata = factor(strata)

  if (is.null(treated.weights)) {
    treated.weights = rep(1 / N1, N1)
  } else {
    stopifnot(length(treated.weights) == N1, all(treated.weights >= 0), sum(treated.weights) > 0)
    treated.weights = treated.weights / sum(treated.weights)
  }

  strata.control = strata[seq_len(N0)]
  strata.treated = strata[(N0 + 1):N]
  treated.levels = unique(as.character(strata.treated))

  fits = list(); rows = list(); dropped = character(0)
  for (s in treated.levels) {
    trt.local = which(as.character(strata.treated) == s)
    ctl.idx   = which(as.character(strata.control) == s)
    share     = sum(treated.weights[trt.local])
    if (length(ctl.idx) < min.controls || share <= 0) {
      if (drop.infeasible) { dropped = c(dropped, s); next }
      stop(sprintf("stratum '%s' has %d treated unit(s) but only %d control(s); need >= %d controls (or set drop.infeasible = TRUE)",
                   s, length(trt.local), length(ctl.idx), min.controls))
    }
    all.idx = c(ctl.idx, N0 + trt.local)
    tw.s    = treated.weights[trt.local] / sum(treated.weights[trt.local])
    fit = synthdid_estimate_weighted(Y[all.idx, , drop = FALSE], length(ctl.idx), T0,
                                     treated.weights = tw.s,
                                     cluster = if (is.null(cluster)) NULL else cluster[all.idx],
                                     ...)
    fits[[s]] = fit
    rows[[s]] = data.frame(stratum = s, N1 = length(trt.local), N0 = length(ctl.idx),
                           share = share, estimate = as.numeric(fit),
                           n1_eff = 1 / sum(tw.s^2),
                           stringsAsFactors = FALSE)
  }
  if (length(fits) == 0) stop("no estimable stratum")

  strata.table = do.call(rbind, rows)
  rownames(strata.table) = NULL
  # Treated-weight mass in dropped strata. Renormalizing the remaining shares
  # changes the target when this is > 0, so callers (and the bootstrap) should
  # record it rather than treat the aggregate as an estimate of the full target.
  dropped.mass = 1 - sum(strata.table$share)
  if (abs(dropped.mass) < 1e-10) dropped.mass = 0   # floating-point residue, not dropped mass
  strata.table$share = strata.table$share / sum(strata.table$share)  # renormalize if any dropped
  estimate = sum(strata.table$share * strata.table$estimate)

  class(estimate) = 'synthdid_estimate_stratified'
  attr(estimate, 'estimator')       = 'synthdid_estimate_stratified'
  attr(estimate, 'strata.fits')     = fits
  attr(estimate, 'strata.table')    = strata.table
  attr(estimate, 'strata.dropped')  = dropped
  attr(estimate, 'dropped.mass')    = dropped.mass
  attr(estimate, 'setup')           = list(Y = Y, N0 = N0, T0 = T0, strata = strata)
  attr(estimate, 'treated.weights') = treated.weights
  attr(estimate, 'cluster')         = cluster
  attr(estimate, 'opts')            = c(list(min.controls = min.controls), list(...))
  estimate
}

#' @method print synthdid_estimate_stratified
#' @export
print.synthdid_estimate_stratified = function(x, ...) {
  cat(sprintf("stratified weighted SDID: %1.5f\n", as.numeric(x)))
  print(attr(x, 'strata.table'))
  invisible(x)
}

#' Calculate Variance-Covariance Matrix for a Stratified Weighted SDID Estimate
#'
#' Bootstrap only. Each replication resamples units (or whole clusters when a
#' cluster vector is available) from the full panel, rebuilds the strata from
#' the resampled rows, and re-runs the entire stratified estimator -- so the
#' within-stratum omega/lambda weights, any detrending, and the stratum shares
#' are all re-fit inside every draw. Strata that lose all treated units or drop
#' below `min.controls` controls in a draw are dropped and the remaining shares
#' renormalized (drop.infeasible = TRUE).
#'
#' Unlike vcov.synthdid_estimate_weighted, regularization parameters are
#' recomputed per draw from the resampled data (full refit) rather than frozen
#' at their original-fit values, since a stratified fit has no single zeta.
#'
#' @param object A synthdid_estimate_stratified model
#' @param method Only "bootstrap" is supported.
#' @param replications the number of bootstrap replications
#' @param cluster An optional vector of cluster IDs; defaults to the cluster stored in the
#'        estimate. When non-NULL, resamples clusters rather than units. Pass cluster = NULL
#'        to force unit-level resampling.
#' @param ... Additional arguments (currently ignored).
#' @method vcov synthdid_estimate_stratified
#' @export
vcov.synthdid_estimate_stratified = function(object,
  method = c("bootstrap"),
  replications = 200,
  cluster = attr(object, 'cluster'), ...) {
    method = match.arg(method)
    samples = stratified_bootstrap_sample(object, replications, cluster)
    samples = samples[!is.na(samples)]
    if (length(samples) < 2) {
      warning("Stratified bootstrap: fewer than 2 valid replicates; returning NA")
      return(matrix(NA_real_))
    }
    R = length(samples)
    matrix((sqrt((R - 1) / R) * sd(samples))^2)
}

# One stratified bootstrap replication for a set of drawn row indices is just a
# full refit; this generates `replications` draws (unit- or cluster-resampled).
# Exposed (not exported) so parallel drivers can call it with their own seeds.
stratified_bootstrap_sample = function(object, replications, cluster = attr(object, 'cluster')) {
  setup = attr(object, 'setup')
  draws = replicate(replications,
    stratified_boot_rep(object, if (is.null(cluster)) "unit" else "cluster", cluster))
  draws
}

#' One stratified bootstrap replication (resample and full refit).
#'
#' Draws a unit- or cluster-level resample of the panel stored in a
#' [synthdid_estimate_stratified] object, rebuilds the strata from the drawn
#' rows, and re-runs the entire stratified estimator with
#' `drop.infeasible = TRUE`. This is the building block of
#' `vcov.synthdid_estimate_stratified`; it is exported so that parallel drivers
#' can call it with their own seeds and reproduce the package's bootstrap exactly.
#'
#' @param object a `synthdid_estimate_stratified` object.
#' @param method `"unit"` or `"cluster"`.
#' @param cluster a vector of cluster IDs of length `nrow(Y)` (required when
#'        `method = "cluster"`).
#' @param return what to return: `"estimate"` (default) gives the refit aggregate
#'        as a numeric scalar; `"detail"` gives `c(estimate, dropped.mass)`, where
#'        `dropped.mass` is the treated-weight mass of strata dropped inside the
#'        draw (so callers can audit how far the resampled target drifts from the
#'        original); `"object"` gives the refit `synthdid_estimate_stratified`
#'        object itself (e.g. to compute a resampled event-study curve).
#' @return see `return`. A draw that empties one side of the panel, or whose
#'         refit fails, yields `NA_real_`, `c(NA, NA)`, or `NULL` respectively.
#' @keywords internal
#' @export
stratified_boot_rep = function(object, method, cluster = NULL,
                               return = c("estimate", "detail", "object")) {
  return = match.arg(return)
  fail = switch(return, estimate = NA_real_,
                detail = c(estimate = NA_real_, dropped.mass = NA_real_),
                object = NULL)
  setup = attr(object, 'setup')
  tw    = attr(object, 'treated.weights')
  opts  = attr(object, 'opts')
  N0 = setup$N0; N = nrow(setup$Y)

  if (method == "cluster") {
    stopifnot(!is.null(cluster), length(cluster) == N)
    cc = cluster[1:N0]; ct = cluster[(N0 + 1):N]
    drawn = sample(unique(cluster), replace = TRUE)
    control.ind   = unlist(lapply(drawn, function(g) which(cc == g)))
    treated.local = unlist(lapply(drawn, function(g) which(ct == g)))
  } else {
    ind = sample(1:N, replace = TRUE)
    control.ind   = sort(ind[ind <= N0])
    treated.local = sort(ind[ind > N0]) - N0
  }
  if (length(control.ind) == 0 || length(treated.local) == 0) return(fail)

  all.ind = c(control.ind, N0 + treated.local)
  args = c(list(Y = setup$Y[all.ind, , drop = FALSE],
                N0 = length(control.ind), T0 = setup$T0,
                strata = setup$strata[all.ind],
                treated.weights = sum_normalize(tw[treated.local]),
                cluster = NULL, drop.infeasible = TRUE),
           opts)
  fit = tryCatch(do.call(synthdid_estimate_stratified, args), error = function(e) NULL)
  if (is.null(fit)) return(fail)
  switch(return,
         estimate = as.numeric(fit),
         detail   = c(estimate = as.numeric(fit), dropped.mass = attr(fit, 'dropped.mass')),
         object   = fit)
}
