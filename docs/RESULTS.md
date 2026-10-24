# Figures 6, 7, and 12

[Overview](../README.md) | [Method](METHOD.md) | [API](CODE.md) | [Validation](REPRODUCIBILITY.md)

Numbering follows [arXiv:2310.03078v2](https://arxiv.org/abs/2310.03078v2), by
Kumar, Aziz, Chakraborty, Ludwig, Gopalakrishnan, Pixley, and Vasseur. These
dual-Haar calculations are the contribution covered by this repository.
The images below are published results, not outputs of the small examples.

## Figure 6: Open-Boundary Free Energy

![Figure 6](figures/fig06.png)

Open-boundary record probabilities provide the finite-size free-energy observable.

[Original PDF](figures/fig06.pdf) | [Research kernel](../research/open_boundary_haar.jl)
| [Runnable example](../examples/quickstart.jl).

The example exercises gate evolution and Born measurements and saves cumulative
record surprisal. Its eight trajectories per size are not a central-charge fit.

## Figure 7: Boundary Dissipation

![Figure 7](figures/fig07.png)

Comparing dissipative and open boundaries probes the boundary contribution.

[Statevector and density-matrix methods](../src/BoundaryDephasing.jl)
| [Boundary-dephasing driver](../examples/dephasing.jl)
| [Algorithm](METHOD.md).

Use the density reference and hidden-history sampling checks before increasing
system size or attempting a boundary-dimension fit. Production-quality estimates
require long-time and ensemble convergence beyond the demonstration settings.

## Figure 12: Monte Carlo Convergence

![Figure 12](figures/fig12.png)

Published Monte Carlo extrapolation and density-matrix comparison.

[Original PDF](figures/fig12.pdf) | [Validation configuration](../examples/validation.toml)
| [Saved numerical checks](../examples/dephasing_expected/).

## Executable Validation

![Actual short-run validation](../examples/dephasing_validation.png)

This separate plot comes from the included code. Both calculations use the same
records; the displayed differences isolate finite hidden-history sampling error.
It tests the method on short circuits and does not recreate the published ensemble.

All paper assets have [source attribution](figures/README.md) and
[checksums](figure_manifest.json).
