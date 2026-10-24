# Boundary Transfer Matrix: Dual-Haar Circuits

[![Tests](https://github.com/ka784-jpg/Boundary-Transfer-Matrix/actions/workflows/tests.yml/badge.svg)](https://github.com/ka784-jpg/Boundary-Transfer-Matrix/actions/workflows/tests.yml)
[![Code license: MIT](https://img.shields.io/badge/code_license-MIT-blue.svg)](LICENSE)

Julia simulations of measurement-record probabilities in monitored dual-unitary
circuits with open and dissipative boundaries. The repository provides statevector
Monte Carlo, an independent full density-matrix reference, finite-sampling
extrapolation, tests, and reproducible examples.

**Paper:** Abhishek Kumar, Kemal Aziz, Ahana Chakraborty, Andreas W. W. Ludwig,
Sarang Gopalakrishnan, J. H. Pixley, and Romain Vasseur,
[*Boundary transfer matrix spectrum of measurement-induced transitions*](https://arxiv.org/abs/2310.03078v2),
Phys. Rev. B **109**, 014303 (2024).
[Published article](https://doi.org/10.1103/PhysRevB.109.014303).

[Run an example](#quickstart) | [Paper figures](docs/RESULTS.md) | [Method](docs/METHOD.md) | [API](docs/CODE.md) | [Validation](docs/REPRODUCIBILITY.md)

![Paper Figure 6: dual-Haar free-energy scaling with open boundaries](docs/figures/fig06.png)

**What the calculation probes.** Measurement-record likelihoods define a
free-energy-like observable. Comparing open and dissipative boundaries isolates
boundary contributions. The gallery covers the dual-Haar contribution to
Figures 6, 7, and 12; other coauthors' circuit families are outside this repository's scope.

## Quickstart

Install [Julia 1.12](https://julialang.org/downloads/), clone the repository, and
run from its root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --startup-file=no examples/dephasing.jl examples/dephasing.toml results/dephasing
julia --project=. --startup-file=no tests/runtests.jl
```

The example samples observed measurement records, estimates their probabilities
by averaging hidden boundary phase-flip histories, and checks them against full
density matrices. It runs at L=4 with both Haar and product initial states.
Outputs include CSV data, error diagnostics, seeds, configuration, runtime version,
and source hashes. Existing output directories are refused. Initial compilation
can take several minutes.

To plot the output with Python 3.12 or newer:

```sh
python -m pip install -r examples/requirements.txt
python examples/plot_dephasing.py results/dephasing results/dephasing/convergence.png
```

### Density-Matrix Validation

![Monte Carlo errors relative to full density matrices](examples/dephasing_validation.png)

This is an actual short-run benchmark, not a published figure: L=2,4,6,8,
12 layers, eight independent records per initial-state ensemble, and up to
4,096 hidden histories per record. Error bars are standard errors of paired
MC-minus-density differences across records.
[Configuration](examples/validation.toml) | [Saved data](examples/dephasing_expected/) | [Interpretation and limits](docs/REPRODUCIBILITY.md).

### Open-Boundary Example

```sh
julia --project=. --startup-file=no examples/quickstart.jl results/open
python examples/plot.py results/open results/open/trace.png
```

![Small-system open-boundary trajectories](examples/trace.png)

This example uses the preserved open-boundary research kernel through a bounded,
seeded adapter. [Saved output](examples/expected/) and [API details](docs/CODE.md).

## Published Results

| Figure 6 | Figure 7 | Figure 12 |
| --- | --- | --- |
| [![Open boundaries](docs/figures/fig06.png)](docs/RESULTS.md#figure-6-open-boundary-free-energy) | [![Dissipative boundary](docs/figures/fig07.png)](docs/RESULTS.md#figure-7-boundary-dissipation) | [![Monte Carlo extrapolation](docs/figures/fig12.png)](docs/RESULTS.md#figure-12-monte-carlo-convergence) |

## Reliability and Scope

- Tests cover dual unitarity, independent dense gate application, Born sampling,
  exact hidden-history sums, density-matrix agreement, log-domain stability,
  reproducible seeds, fits, and output protection.
- Statevector memory grows as 2^L; density-matrix memory grows as 4^L.
  The density reference is capped at L=10; the statevector API at L=20.
- The examples validate the method on small systems. They do not regenerate the
  paper's large production ensembles or certify its fitted critical parameters.
- [Method documentation](docs/METHOD.md) explains why probabilities must be
  averaged **before** taking their logarithm, and how to inspect sampling error.

Code is [MIT licensed](LICENSE). Paper figures are excluded from MIT and retain
their original rights and [attribution](docs/figures/README.md). Cite the full paper.
