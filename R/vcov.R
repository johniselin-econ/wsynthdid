#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates based on the following three options
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item The jackknife, Algorithm 3 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#' }
#'
#' The jackknife is not recommended for SC, see section 5 in Arkhangelsky et al.
#' "placebo" is the only option that works for only one treated unit.
#'
#' @param object A synthdid model
#' @param method, the CI method. The default is bootstrap (warning: this may be slow on large
#'  data sets, the jackknife option is the fastest, with the caveat that it is not recommended
#'  for SC).
#' @param replications, the number of bootstrap replications
#' @param ... Additional arguments (currently ignored).
#'
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
#'  "Synthetic Difference in Differences". arXiv preprint arXiv:1812.09970, 2019.
#'
#' @method vcov synthdid_estimate
#' @export
vcov.synthdid_estimate = function(object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200, ...) {
    method = match.arg(method)
    if(method == 'bootstrap') {
	se = bootstrap_se(object, replications)
    } else if(method == 'jackknife') {
	se = jackknife_se(object)
    } else if(method == 'placebo') {
	se = placebo_se(object, replications)
    }
    matrix(se^2)
}

#' Calculate the standard error of a synthetic diff in diff estimate. Deprecated. Use vcov.synthdid_estimate.
#' @param ... Any valid arguments for vcov.synthdid_estimate
#' @export synthdid_se
synthdid_se = function(...) { sqrt(vcov(...)) }


# The bootstrap se: Algorithm 2 of Arkhangelsky et al.
bootstrap_se = function(estimate, replications) { sqrt((replications-1)/replications) * sd(bootstrap_sample(estimate, replications)) }
bootstrap_sample = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
    theta = function(ind) {
	if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
	else {
	    weights.boot = weights
	    weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
	    do.call(synthdid_estimate, c(list(Y=setup$Y[sort(ind),], N0=sum(ind <= setup$N0), T0=setup$T0, X=setup$X[sort(ind), ,], weights=weights.boot), opts))
	}
    }
    bootstrap.estimates = rep(NA, replications)
    count = 0
    while(count < replications) {
	bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
	if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
    }
    bootstrap.estimates
}


# The fixed-weights jackknife estimate of variance: Algorithm 3 of Arkhangelsky et al.
# if weights = NULL is passed explicitly, calculates the usual jackknife estimate of variance.
# returns NA if there is one treated unit or, for the fixed-weights jackknife, one control with nonzero weight
jackknife_se = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if (!is.null(weights)) {
      opts$update.omega = opts$update.lambda = FALSE
    }
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) { return(NA) }
    theta = function(ind) {
	weights.jk = weights
	if (!is.null(weights)) { weights.jk$omega = sum_normalize(weights$omega[ind[ind <= setup$N0]]) }
	estimate.jk = do.call(synthdid_estimate,
	    c(list(Y=setup$Y[ind, ], N0=sum(ind <= setup$N0), T0=setup$T0, X = setup$X[ind, , ], weights = weights.jk), opts))
    }
    jackknife(1:nrow(setup$Y), theta)
}

#' Jackknife standard error of function `theta` at samples `x`.
#' @param x vector of samples
#' @param theta a function which returns a scalar estimate
#' @importFrom stats var
#' @keywords internal
jackknife = function(x, theta) {
  n = length(x)
  u = rep(0, n)
  for (i in 1:n) {
    u[i] = theta(x[-i])
  }
  jack.se = sqrt(((n - 1) / n) * (n - 1) * var(u))

  jack.se
}



# The placebo se: Algorithm 4 of Arkhangelsky et al.
placebo_se = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    N1 = nrow(setup$Y) - setup$N0
    if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
    theta = function(ind) {
	N0 = length(ind)-N1
	weights.boot = weights
	weights.boot$omega = sum_normalize(weights$omega[ind[1:N0]])
        do.call(synthdid_estimate, c(list(Y=setup$Y[ind,], N0=N0,  T0=setup$T0,  X=setup$X[ind, ,], weights=weights.boot), opts))
    }
    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}

#' Normalize a non-negative weight vector to sum to one.
#'
#' Returns `x / sum(x)`, or uniform weights when `x` sums to zero. Used for the
#' weight renormalization inside every resampling variance estimator, and
#' exported so that external parallel drivers (e.g. a replication script's own
#' bootstrap loop) can reproduce the package's renormalization exactly.
#'
#' @param x a numeric vector of non-negative weights.
#' @return a numeric vector of the same length summing to one.
#' @keywords internal
#' @export
sum_normalize = function(x) {
    if(sum(x) != 0) { x / sum(x) }
    else { rep(1/length(x), length(x)) }
    # if given a vector of zeros, return uniform weights
    # this fine when used in bootstrap and placebo standard errors, where it is used only for initialization
    # for jackknife standard errors, where it isn't, we handle the case of a vector of zeros without calling this function.
}
# =============================================================================
# WEIGHTED VERSIONS OF VARIANCE ESTIMATION
# =============================================================================

#' Calculate Variance-Covariance Matrix for a Weighted Synthetic DID Estimate
#'
#' Provides variance estimates for weighted synthdid estimates.
#' The key difference from the unweighted version is that treated unit weights
#' must be renormalized when resampling.
#'
#' @param object A synthdid_estimate_weighted model
#' @param method, the CI method. The default is bootstrap.
#' @param replications, the number of bootstrap replications
#' @param placebo.weights For the placebo method, how to weight pseudo-treated units:
#'        "uniform" (default, 1/N1 for each), "size_match" (weight by pre-treatment outcome levels),
#'        or "permute" (randomly assign the original treated.weights vector).
#' @param cluster An optional vector of cluster IDs for cluster-robust inference. If NULL (default),
#'        uses the cluster stored in the estimate object (if any). When non-NULL, resamples clusters
#'        rather than individual units.
#' @param ... Additional arguments (currently ignored).
#'
#' @method vcov synthdid_estimate_weighted
#' @export
vcov.synthdid_estimate_weighted = function(object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200,
  placebo.weights = c("uniform", "size_match", "permute"),
  cluster = attr(object, 'cluster'), ...) {
    method = match.arg(method)
    placebo.weights = match.arg(placebo.weights)

    if (!is.null(cluster)) {
      # Cluster-robust inference
      if(method == 'bootstrap') {
        se = cluster_bootstrap_se_weighted(object, replications, cluster)
      } else if(method == 'jackknife') {
        se = cluster_jackknife_se_weighted(object, cluster)
      } else if(method == 'placebo') {
        se = placebo_se_weighted(object, replications, placebo.weights)
        warning("Placebo SE does not currently support clustering; using unit-level placebo.")
      }
    } else {
      # Unit-level inference (existing behavior)
      if(method == 'bootstrap') {
        se = bootstrap_se_weighted(object, replications)
      } else if(method == 'jackknife') {
        se = jackknife_se_weighted(object)
      } else if(method == 'placebo') {
        se = placebo_se_weighted(object, replications, placebo.weights)
      }
    }
    matrix(se^2)
}

#' Calculate the standard error of a weighted synthetic diff in diff estimate
#' @param ... Any valid arguments for vcov.synthdid_estimate_weighted
#' @export synthdid_se_weighted
synthdid_se_weighted = function(...) { sqrt(vcov(...)) }


# The bootstrap se for weighted estimates: modified Algorithm 2
# Key change: renormalize treated.weights for resampled treated units
bootstrap_se_weighted = function(estimate, replications) {
  samples = bootstrap_sample_weighted(estimate, replications)
  if (length(samples) == 0 || all(is.na(samples))) {
    warning("Weighted bootstrap: no valid replicates; returning NA for SE")
    return(NA)
  }
  sqrt((replications-1)/replications) * sd(samples, na.rm = TRUE)
}

bootstrap_sample_weighted = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    treated.weights = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N1 = nrow(setup$Y) - setup$N0

    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }

    theta = function(ind) {
      control.ind = sort(ind[ind <= setup$N0])
      treated.ind = sort(ind[ind > setup$N0])
      treated.ind.local = treated.ind - setup$N0  # indices within treated units

      if(length(control.ind) == 0 || length(treated.ind) == 0) { return(NA) }

      # Renormalize control weights
      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[control.ind])

      # Assign weights to each resampled treated row (one weight per row, matching Y.boot)
      treated.weights.boot = treated.weights[treated.ind.local]
      treated.weights.boot = sum_normalize(treated.weights.boot)

      # Reconstruct Y matrix with resampled units
      Y.boot = setup$Y[c(control.ind, treated.ind), ]
      X.boot = setup$X[c(control.ind, treated.ind), , ]
      N0.boot = length(control.ind)

      do.call(synthdid_estimate_weighted,
              c(list(Y = Y.boot, N0 = N0.boot, T0 = setup$T0, X = X.boot,
                     treated.weights = treated.weights.boot,
                     period.weights = period.weights,
                     weights = weights.boot), opts))
    }

    bootstrap.estimates = rep(NA, replications)
    count = 0
    max_attempts = replications * 10
    attempts = 0
    while(count < replications && attempts < max_attempts) {
      attempts = attempts + 1
      bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
      if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
    }
    if (count < replications) {
      warning(sprintf("Weighted bootstrap: only %d of %d replicates completed after %d attempts",
                       count, replications, attempts))
    }
    bootstrap.estimates[1:count]
}


# The fixed-weights jackknife estimate for weighted estimates: modified Algorithm 3
# Key change: renormalize treated.weights when leaving out treated units
jackknife_se_weighted = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    treated.weights.orig = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N1 = nrow(setup$Y) - setup$N0

    if (!is.null(weights)) {
      opts$update.omega = opts$update.lambda = FALSE
    }
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) {
      return(NA)
    }

    theta = function(ind) {
      control.ind = ind[ind <= setup$N0]
      treated.ind = ind[ind > setup$N0]
      treated.ind.local = treated.ind - setup$N0

      # Renormalize control weights
      weights.jk = weights
      if (!is.null(weights)) {
        weights.jk$omega = sum_normalize(weights$omega[control.ind])
      }

      # Renormalize treated weights for remaining treated units
      treated.weights.jk = treated.weights.orig[treated.ind.local]
      treated.weights.jk = sum_normalize(treated.weights.jk)

      estimate.jk = do.call(synthdid_estimate_weighted,
          c(list(Y = setup$Y[ind, ], N0 = sum(ind <= setup$N0), T0 = setup$T0,
                 X = setup$X[ind, , ],
                 treated.weights = treated.weights.jk,
                 period.weights = period.weights,
                 weights = weights.jk), opts))
    }
    jackknife(1:nrow(setup$Y), theta)
}


# The placebo se for weighted estimates: modified Algorithm 4
# For placebo, we reassign N1 control units as "treated" and estimate effect
# Key change: need to assign weights to the placebo treated units
# placebo.weights controls how pseudo-treated units are weighted:
#   "uniform"    - equal weights 1/N1 (default, tests sharp null under uniform weighting)
#   "size_match" - weight by pre-treatment outcome levels (mimics size-based weighting)
#   "permute"    - randomly permute the original treated.weights (preserves weight concentration)
placebo_se_weighted = function(estimate, replications, placebo.weights = "uniform") {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    treated.weights.orig = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N1 = nrow(setup$Y) - setup$N0

    if (setup$N0 <= N1) {
      stop('must have more controls than treated units to use the placebo se')
    }

    theta = function(ind) {
      N0.placebo = length(ind) - N1
      placebo.treated.ind = ind[(N0.placebo + 1):length(ind)]

      # Renormalize control weights for remaining controls
      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[ind[1:N0.placebo]])

      # Determine placebo treated weights based on method
      if (placebo.weights == "uniform") {
        tw.placebo = rep(1 / N1, N1)
      } else if (placebo.weights == "permute") {
        tw.placebo = treated.weights.orig[sample.int(N1)]
      } else if (placebo.weights == "size_match") {
        # Weight pseudo-treated units by their pre-treatment outcome levels
        Y.pre = setup$Y[placebo.treated.ind, 1:setup$T0, drop = FALSE]
        w = abs(rowMeans(Y.pre))
        if (sum(w) == 0) w = rep(1, N1)
        tw.placebo = w / sum(w)
      }

      do.call(synthdid_estimate_weighted,
              c(list(Y = setup$Y[ind, ], N0 = N0.placebo, T0 = setup$T0,
                     X = setup$X[ind, , ],
                     treated.weights = tw.placebo,
                     period.weights = period.weights,
                     weights = weights.boot), opts))
    }

    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}


# =============================================================================
# CLUSTER-ROBUST VARIANCE ESTIMATION
# Cluster bootstrap in the spirit of Clarke et al. (2023) / the Stata sdid
# package (Daniel-Pailanir). NOTE: this implementation REFITS omega and lambda
# in every draw (opts keep update.omega = update.lambda = TRUE; the subset
# weights are only the Frank-Wolfe starting point). The cluster jackknife below
# is the fixed-weight variant.
# =============================================================================

# Cluster bootstrap SE: resample clusters with replacement, renormalize the
# subset weights as the starting point, refit
cluster_bootstrap_se_weighted = function(estimate, replications, cluster) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    treated.weights = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N0 = setup$N0
    N = nrow(setup$Y)

    cluster_control = cluster[1:N0]
    cluster_treated = cluster[(N0+1):N]
    unique_clusters = unique(cluster)

    theta = function() {
      drawn = sample(unique_clusters, replace = TRUE)
      control.ind = unlist(lapply(drawn, function(cl) which(cluster_control == cl)))
      treated.ind.local = unlist(lapply(drawn, function(cl) which(cluster_treated == cl)))
      if (length(control.ind) == 0 || length(treated.ind.local) == 0) return(NA)

      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[control.ind])
      tw.boot = sum_normalize(treated.weights[treated.ind.local])

      all.ind = c(control.ind, N0 + treated.ind.local)
      Y.boot = setup$Y[all.ind, , drop = FALSE]
      X.boot = setup$X[all.ind, , , drop = FALSE]
      N0.boot = length(control.ind)

      c(do.call(synthdid_estimate_weighted,
                c(list(Y = Y.boot, N0 = N0.boot, T0 = setup$T0, X = X.boot,
                       treated.weights = tw.boot, period.weights = period.weights,
                       weights = weights.boot), opts)))
    }

    estimates = rep(NA, replications)
    count = 0; max_attempts = replications * 10; attempts = 0; failures = 0
    while (count < replications && attempts < max_attempts) {
      attempts = attempts + 1
      tryCatch({
        est = theta()
        if (!is.na(est)) { count = count + 1; estimates[count] = est }
      }, error = function(e) { failures <<- failures + 1 })
    }
    if (failures > 0) warning(sprintf("Cluster bootstrap: %d of %d attempts failed", failures, attempts))
    if (count < replications) warning(sprintf("Cluster bootstrap: only %d of %d replicates completed", count, replications))
    sqrt((replications-1)/replications) * sd(estimates[1:count])
}

# Cluster jackknife SE: leave out one cluster at a time
cluster_jackknife_se_weighted = function(estimate, cluster) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    treated.weights.orig = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N0 = setup$N0
    N = nrow(setup$Y)

    opts$update.omega = FALSE
    opts$update.lambda = FALSE

    cluster_control = cluster[1:N0]
    cluster_treated = cluster[(N0+1):N]
    unique_clusters = unique(cluster)
    K = length(unique_clusters)
    if (K <= 1) return(NA)

    theta_k = rep(NA, K)
    for (k in 1:K) {
      cl = unique_clusters[k]
      keep_control = which(cluster_control != cl)
      keep_treated = which(cluster_treated != cl)
      if (length(keep_control) == 0 || length(keep_treated) == 0) { next }

      weights.jk = weights
      weights.jk$omega = sum_normalize(weights$omega[keep_control])
      tw.jk = sum_normalize(treated.weights.orig[keep_treated])

      keep.ind = c(keep_control, N0 + keep_treated)
      theta_k[k] = tryCatch(
        c(do.call(synthdid_estimate_weighted,
                  c(list(Y = setup$Y[keep.ind, , drop = FALSE],
                         N0 = length(keep_control), T0 = setup$T0,
                         X = setup$X[keep.ind, , , drop = FALSE],
                         treated.weights = tw.jk, period.weights = period.weights,
                         weights = weights.jk), opts))),
        error = function(e) NA)
    }

    valid = !is.na(theta_k)
    if (sum(valid) < 2) return(NA)
    K.valid = sum(valid)
    sqrt(((K.valid - 1) / K.valid) * (K.valid - 1) * var(theta_k[valid]))
}
