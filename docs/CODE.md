# Code Interfaces

```julia
include("examples/kernels.jl")
surprisal = OpenBoundary.trajectory(; L=4, p=0.14, steps=64, seed=2026)
```

The result is cumulative measurement-record surprisal at each layer, in nats.
Its ensemble expectation is the Shannon entropy of the record. A trajectory's
value alone is not a free-energy estimate with controlled statistical error.

The wrapper calls the archived `rangate`, `gate`, `even_projmsm`, and
`odd_projmsm` definitions. It preserves the open-boundary schedule: the first
alternating gate layer pairs (1,2), (3,4), etc.; the other pairs (2,3), (4,5), etc.
Measurements act on all sites after the first kind of layer and on interior
sites after the other kind, as in the preserved template. There is no wrapping
bond. The state is initially Haar-random, represented by a normalized complex
Gaussian vector.

The example validates even L in 2:10, finite p in [0,1], layer count in 1:320,
and a nonnegative 32-bit seed. Calls reset module globals and the default RNG.
It is not a thread-safe API; independent workers need separate processes.

## Archive Versus Adapter

[research/open_boundary_haar.jl](../research/open_boundary_haar.jl) is the
unaltered research template. It retains a top-level batch launch. The adapter
uses Julia's parsed-expression `include` hook to omit that launch and original
import declarations, supplying the runtime imports explicitly. In particular,
`Normal` is supplied by Distributions; the original file did not import it
explicitly. The local gate and measurement function bodies are unchanged.
Unused CSV, DataFrames, and GLM imports are unnecessary for this example.

Do not run archived sources as the quickstart. The supported commands are in
the [README](../README.md). The archive contains old helper functions that the
current tests do not certify, including an unused entropy helper with a legacy
global-name dependency.

The original boundary-dephasing Monte Carlo and density-matrix codes for
Figures 7 and 12 have not yet been recovered. The presence of those published
images does not imply that this adapter implements them.
