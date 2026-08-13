# Example scripts

Sixteen self-contained R scripts that exercise the whole public API of `neuralsbi`. Each one runs on its own, from a fresh session, in a few minutes.

Every simulator is taken from a published example in another package or tutorial, so the model is somebody else's and the inference is ours. The source, with a link, is in the header of each script.

| script | what it covers | model source | runtime |
| --- | --- | --- | --- |
| `00_installation.R` | installing the package and the torch back end, optional dependencies, a smoke test | none | 5 s |
| `01_basic_npe_example.R` | the shortest end-to-end NPE run on real data | SIR fit to the 1978 boarding-school influenza outbreak, [Stan case study](https://mc-stan.org/learn-stan/case-studies/boarding_school_case_study.html) | 2 min |
| `02_priors.R` | `prior_uniform()`, `prior_normal()`, `prior_custom()`, `sample_prior()`, `within_support()` | Sick-Sicker calibration ranges, [darthpack](https://github.com/DARTH-git/darthpack) | 1 s |
| `03_npe.R` | the four density estimators, training controls, `save_npe()`, `map_estimate()`, `log_prob()`, `npe_sequential()` | Sick-Sicker cohort model, [DARTH cSTM tutorial](https://github.com/DARTH-git/cohort-modeling-tutorial-intro) | 5 min |
| `04_nle.R` | `nle()`, `log_lik()`, `likelihood_fn()`, MCMC posteriors, conditioning on n i.i.d. observations | one-compartment PK model for `Theoph`, Pinheiro and Bates via `nlme::SSfol` | 4 min |
| `05_nre.R` | `nre()`, `log_ratio()`, the four classifiers, `num_atoms` | stochastic SIR, [mcstate SIR vignette](https://mrc-ide.github.io/mcstate/articles/sir_models.html) | 3 min |
| `06_embedding_networks.R` | `embedding_mlp()`, choosing `output_dim` | SIR fit to the 1948 Consett measles outbreak, [SBIED lesson 2](https://kingaa.github.io/sbied/stochsim/) | 5 min |
| `07_npe_diagnostics.R` | the whole diagnostic suite on one fit, and what each failure means | boarding-school influenza SIR, as in `01` | 5 min |
| `08_nle_with_stan.R` | `stan_code()`, `stan_data()`, `write_stan_model()`, `sampler = "stan"` | BCG vaccine meta-analysis, [metafor / metadat](https://wviechtb.github.io/metadat/reference/dat.colditz1994.html) | 3 min |
| `09_sbc.R` | `sbc()`, `expected_coverage()`, what a failing rank histogram looks like, and the failure SBC cannot see | `task_gaussian_linear()` and `task_two_moons()`, sbibm | 2 min |
| `10_c2st.R` | `c2st()`, what the number means and what it cannot see | `task_gaussian_linear()`, sbibm | 5 min |
| `11_tarp.R` | `tarp()`, joint versus marginal calibration, reference points | `task_gaussian_linear()` and `task_two_moons()`, sbibm | 1 min |
| `12_plots.R` | `pairplot()`, `plot_sbc()`, `plot_coverage()`, `plot_tarp()`, `plot_posterior_predictive()` | `task_sir()`, sbibm | 5 min |
| `13_prior_posterior_predictive_checks.R` | prior and posterior predictive checks, and the failure a predictive check cannot catch | HIV therapy Markov model, [heemod vignette](https://cran.r-project.org/web/packages/heemod/vignettes/c_homogeneous.html) | 3 min |
| `14_utilities.R` | the simulator contract, `sim_args`, `future` plans, progress, seeds, saving, tidy accessors | `task_gaussian_linear()` and the Sick-Sicker model | 35 s |
| `15_npe_vs_pomp.R` | NPE against iterated filtering on the same POMP | Consett measles SIR, [SBIED lesson 2](https://kingaa.github.io/sbied/stochsim/) | 3 min |

Runtimes are wall clock on four CPU cores with libtorch installed.

## Running them

```sh
Rscript dev/01_basic_npe_example.R
```

Scripts that need `torch` say so and fall back to `density_estimator = "linear_gaussian"` where a fallback makes sense. `06_embedding_networks.R` and `12_plots.R` stop early without `torch` and `ggplot2` respectively, because there is nothing left to demonstrate.

Optional packages used, none of them required by the package itself: `ggplot2`, `GGally` and `ggdensity` for plots, `future` for parallel simulation, `progressr` for progress bars, `cmdstanr` or `rstan` for the Stan half of `08`, and `pomp` for `15`.

## A note on training controls

Several scripts pass `max_epochs` and `patience` below their defaults so they finish in a few minutes. Where that happens the script says so. Use the defaults for real work.

## A note on the log transform

`05`, `06` and `15` hand the estimator `log(1 + counts)` rather than the counts. The reasoning, the measured size of the effect (smaller than one comparison suggests), and what other packages do about it are in the header of `05_nre.R` and in issue #173. `12_plots.R` needs no such transform because `task_sir()` emits binomial proportions rather than unbounded counts, following sbibm.
