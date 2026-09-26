# wsynthdid: Weighted Synthetic Difference-in-Differences

This is a fork of the [`synthdid`](https://github.com/synth-inference/synthdid) R package that extends the synthetic difference-in-differences estimator of Arkhangelsky et al. (2021) with **user-specified weights on treated units**.

The standard SDID estimator averages treated units equally. When treated units differ in size---for example, counties ranging from 10,000 to 10 million residents---the equally-weighted and population-weighted estimands can diverge substantially. This package allows the researcher to specify weights that match their estimand of interest.

### What this fork adds

- **`synthdid_estimate_weighted()`**: weighted SDID with researcher-specified treated-unit weights
- **`sc_estimate_weighted()`**, **`did_estimate_weighted()`**: weighted SC and DiD variants
- **`synthdid_estimate_staggered()`**: staggered-adoption weighted SDID (cohort-wise estimation against not-yet-treated or never-treated controls, aggregated with researcher-chosen cohort weights)
- **`synthdid_estimate_stratified()`**: weighted SDID with the donor pool restricted to researcher-defined strata
- **`synthdid_event_study()`**: event-study decomposition for both weighted and unweighted estimates, with unit-level or cluster-robust bootstrap confidence bands
- **Adapted variance estimators**: bootstrap, jackknife, and placebo SEs with weight renormalization
- **Cluster-robust standard errors**: cluster bootstrap and jackknife for settings where treatment is assigned at a higher level than the unit of observation
- **Placebo weight options**: uniform, size-matched, and permuted weights for placebo inference
- **Effective sample size tuning**: optional Kish's N<sub>eff</sub> adjustment for the regularization parameter

### Installation

```R
devtools::install_github("johniselin-econ/wsynthdid")          # latest
devtools::install_github("johniselin-econ/wsynthdid@v0.1.1")   # the version used in the paper
```

### Quick example

```R
library(wsynthdid)

# Standard SDID (equal weights)
data('california_prop99')
setup = panel.matrices(california_prop99)
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
sprintf('point estimate: %1.2f', tau.hat)

# Weighted SDID (custom weights on treated units)
N1 = nrow(setup$Y) - setup$N0
pop_weights = runif(N1)  # replace with real population weights
tau.w = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
          treated.weights = pop_weights)
sprintf('weighted estimate: %1.2f', tau.w)

# Cluster-robust standard errors
cluster = rep(1:5, length.out = nrow(setup$Y))
tau.cl = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
           treated.weights = pop_weights, cluster = cluster)
se = sqrt(vcov(tau.cl, method = "bootstrap", replications = 200))

# Event study (se.method = "bootstrap" adds CIs; an estimate's stored
# cluster is used automatically, so tau.cl would get cluster-robust bands)
es = synthdid_event_study(tau.w)
plot_event_study(es)
```

### Paper

John Iselin and Erica Ryan. **Weighted Synthetic Difference-in-Differences**. Working paper, 2026.

The paper applies the weighted estimator to the county-level mortality effects of the 2014 ACA Medicaid expansion, following Borgschulte and Vogler (2020). The equally-weighted SDID estimate is near zero. The unconditional population-weighted estimate is about −17.5 deaths per 100,000 but fails an in-time placebo: mortality trends vary with county size, and control weights re-solved against the population-weighted target remain size-unrepresentative. The size-stratified estimator (`synthdid_estimate_stratified()`), which keeps the population-weighted estimand while comparing each county only with same-size donors, gives about −4.8 (SE 3.9) with no detectable all-cause placebo effect; its interval admits both meaningful reductions and zero. The weighting choice selects both the estimand and the identification burden.

The manuscript, supplement, and full replication package (data pipeline, application analysis, Monte Carlo sweeps, and an `renv` lockfile that pins this package at v0.1.1, commit `3c25415`) live in a separate repository, [**johniselin-econ/weighted-sdid**](https://github.com/johniselin-econ/weighted-sdid), which will be made public with the paper. This repository is the estimator only.

### Citation

See [`CITATION.cff`](CITATION.cff), or GitHub's "Cite this repository" button. Please also cite Arkhangelsky et al. (2021) for the underlying estimator.

### References

Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
**Synthetic Difference in Differences.** *American Economic Review*, 111(12):4088-4118, 2021.

Damian Clarke, Daniel Pailanir, Susan Athey, and Guido Imbens.
**Synthetic Difference-in-Differences Estimation.** IZA Discussion Paper No. 15907, 2023.

Mark Borgschulte and Jacob Vogler.
**Did the ACA Medicaid Expansion Save Lives?** *Journal of Health Economics*, 72:102333, 2020.

### Upstream

This fork is based on [synth-inference/synthdid](https://github.com/synth-inference/synthdid). The original package API is fully preserved; all weighted functions are additive extensions.

### Disclosure

The authors used Claude (Anthropic) to assist with code development, data processing, and manuscript editing. All analyses and results were reviewed and verified by the authors.
