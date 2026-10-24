# Simulation Interfaces

[Overview](../README.md) | [Method](METHOD.md) | [Validation](REPRODUCIBILITY.md)

## Boundary-Dephasing Statevectors

```julia
using Random
include("src/BoundaryDephasing.jl")
using .BoundaryDephasing
rng = Xoshiro(2026)
circuit = random_circuit(rng; L=4, p=0.14, layers=12)
psi = initial_state(rng, 4; kind=:haar)
record = sample_record(rng, circuit, psi)
mc = monte_carlo(Xoshiro(99), circuit, psi, record; counts=[64, 256, 1024])
reference = density_record(circuit, psi, record)
f = free_energy_density(mc.logprob[:, end], 4)
```

The observed record is sampled with its Born probability. The independent inner
RNG estimates its marginal likelihood over unobserved boundary jumps.

| Function | Result / contract |
| --- | --- |
| `Layer(bonds, gates, measured; dephase=false)` | Ordered two-site gates, measured sites, optional final boundary dephasing |
| `Circuit(L, layers)` | Validates even L in 2:20, unitary 4x4 gates, and site indices |
| `random_circuit(rng; L,p,layers,dephase)` | Open brickwork circuit; dephasing after each even half-layer |
| `initial_state(rng,L;kind)` | Normalized `:haar` or computational-basis `:product` state |
| `sample_record(rng,circuit,psi)` | One list of observed binary outcomes per layer |
| `path_logprob(circuit,psi,record,jumps)` | Prefix log probabilities conditional on a fixed I/Z history |
| `monte_carlo(rng,circuit,psi,record;counts)` | Prefix log probabilities, relative probability SE, and effective sample size at each cumulative count |
| `density_record(circuit,psi,record)` | Independent full-density prefix likelihood and final conditional density; L <= 10 |
| `enumerate_histories(circuit,psi,record)` | Exact uniform I/Z-history sum; at most 16 dephasing events |
| `free_energy_density(logprob,L;burnin=0)` | Record surprisal per site and retained layer, in nats |
| `extrapolate_samples(counts,values;degree=4,min_samples=300)` | Polynomial intercept at infinite sample count, condition number, residual |
| `bootstrap_extrapolation(rng,counts,curves;replicates=500)` | Resamples independent record curves, preserving correlations between checkpoints |
| `fit_boundary_dimension(sizes,deltaf;Lmin)` | Unweighted fit of L delta-f to a + pi h/L |

One time unit is a gate-plus-measurement half-layer. Site 1 is the least
significant computational-basis bit. Measurements follow gates; boundary
dephasing follows measurements. Invalid configurations fail explicitly.
No clipping is used to conceal a zero Monte Carlo likelihood: it remains
`-Inf`, and free-energy conversion rejects non-finite inputs.

All stochastic APIs accept an explicit RNG. Separate calls with distinct RNGs
and unmodified circuit inputs do not share simulation state. The supplied
ensemble driver is serial and records its seed rules in `metadata.toml`.

## Open-Boundary Research Kernel

```julia
include("examples/kernels.jl")
surprisal = OpenBoundary.trajectory(; L=4, p=0.14, steps=64, seed=2026)
```

[research/open_boundary_haar.jl](../research/open_boundary_haar.jl) is preserved
byte-for-byte, with a [SHA-256 manifest](source_manifest.json). The adapter loads
its gate and measurement definitions without invoking its batch launcher.
It supplies the required runtime imports and skips unused imports.

This adapter uses module globals and resets the default RNG. It is serial,
supports even L in 2:10 and 1:320 layers, and requires separate processes for
independent concurrent workers. Its unused historical helpers are not certified
by the tests. Use the documented driver, not the template's top-level launcher.
