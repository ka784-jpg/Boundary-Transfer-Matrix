# Reproducibility and Recovery Status

The supported example reuses the original open-boundary dynamics and records
its source hash, Julia version, and explicit seed rule. It does not reproduce
the high-statistics production calculation or estimate critical exponents.

The source was recovered on 2026-10-24 from the author's Amarel archive,
`old/ceff_openbc_fulltime_in_f/inputs/spins8_haar`, beneath the home directory.
On that host `/home/ka784` resolves to `/cache/home/ka784`. A descriptive local
filename was assigned without modifying file bytes. See
[source_manifest.json](source_manifest.json).

The user identified Figures 6, 7, and 12 as their contribution and confirmed
that the dissipative calculation used statevectors with Monte Carlo sampling.
The original Fig. 7/12 sources have not yet been identified in the inspected
Amarel/AP41 directories. Original production seeds, full ensembles, final fit
scripts, and a historical environment lock are also unavailable here. Recovery
of one related template does not establish an exact production commit.

The new examples target Julia 1.12.6 with a pinned Manifest. Floating-point
rounding and RNG behavior can change across versions and platforms. Tests use
physical invariants and same-process seeded repeatability, rather than requiring
bitwise equality to saved outputs. The example uses natural logarithms and
defines a time step as one gate-plus-measurement layer.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --startup-file=no tests/runtests.jl
python tests/integrity.py
```

The small dense checks validate local kernels, not the unrecovered Monte Carlo
method. Full Fig. 7/12 reproduction must wait for the original code and data;
an independently written replacement must be labeled as such.
