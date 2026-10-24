# Boundary-Dephasing Method

[API](CODE.md) | [Validation](REPRODUCIBILITY.md) | [Paper, Appendix D](https://arxiv.org/html/2310.03078v2)

## Circuit and Channel

Open brickwork layers alternate bonds (1,2),(3,4),... with (2,3),(4,5),... .
The first type measures all sites independently with probability p; the second
measures interior sites. Complete dephasing acts on the rightmost site after
each second layer:

$$\mathcal D(\rho)=P_0\rho P_0+P_1\rho P_1=\tfrac12(\rho+Z_L\rho Z_L).$$

The gate sampler uses independent Haar U(2) factors around
`exp[-i*pi/4*(XX+YY)-i*J*ZZ]`, with J uniform on [0,2*pi), matching the
open-boundary research kernel's convention. An independent global phase has no
effect on measurement probabilities. Gate and reshuffled-gate unitarity are tested.

## Hidden Histories and Observed Records

For fixed gates, measurement locations, and observed record m, sample each
unobserved dephasing event as I or Z with equal probability. If z is one such
history, accumulate `log P(m|z)` by normalizing after each recorded projection.
The observed likelihood is

$$P(m)=2^{-N_D}\sum_z P(m\mid z),\qquad
\widehat P_\tau(m)=\frac1\tau\sum_{a=1}^{\tau}P(m\mid z_a).$$

I and Z are normalized unitaries. Their sampling probability already supplies
the channel's factor of one half; adding an extra Kraus normalization would
double-count that weight. The code averages likelihoods using log-sum-exp and
only then computes minus their logarithm. Averaging conditional surprisals
instead would answer a different question involving the unobserved environment.

The outer ensemble draws observed records using physical Born measurements and
hidden phase flips. The inner likelihood estimator uses a separate RNG stream.
For one record, its retained-window free-energy contribution is
`-[log P(m_1:T)-log P(m_1:b)]/[L*(T-b)]`, where b is the burn-in layer count.
Ensemble averages of this quantity estimate the corresponding record entropy rate.

## Sampling and Fits

Finite sampling introduces bias after taking a logarithm. The API supports a
fourth-order fit in inverse sample count and whole-record bootstrap resampling.
The example includes checkpoints at 300 through 600 as well as larger counts;
all selected counts >=300 enter its fit. To restrict the window to 300:600,
pass only those columns to the fitting function.

Inspect relative probability errors, effective sample sizes, polynomial condition
numbers, and the stability of the extrapolated intercept. A quartic fit can be
unstable on small ensembles; a fit's existence is not evidence of convergence.
Zero estimated probabilities remain explicit failures rather than being replaced
by arbitrary positive floors. Longer circuits may require substantially more
hidden histories than the examples.

For size scaling, `fit_boundary_dimension` fits
`L*(f_dissipative-f_open)=surface_difference+pi*h/L`. It is an unweighted helper,
not a substitute for uncertainty propagation, equilibration, or fit-window checks.
