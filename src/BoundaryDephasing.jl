module BoundaryDephasing

using LinearAlgebra, Random, SparseArrays, Statistics

export Layer, Circuit, random_circuit, initial_state, sample_record,
       path_logprob, density_record, enumerate_histories, monte_carlo,
       free_energy_density, extrapolate_samples, bootstrap_extrapolation,
       fit_boundary_dimension

struct Layer
    bonds::Vector{Int}
    gates::Vector{Matrix{ComplexF64}}
    measured::Vector{Int}
    dephase::Bool
end
Layer(bonds, gates, measured; dephase=false) =
    Layer(Int.(collect(bonds)), Matrix{ComplexF64}.(gates), Int.(collect(measured)), dephase)

struct Circuit
    L::Int
    layers::Vector{Layer}
    function Circuit(L::Integer, layers::Vector{Layer})
        iseven(L) && 2 <= L <= 20 || throw(ArgumentError("L must be even, from 2 to 20"))
        !isempty(layers) || throw(ArgumentError("at least one layer is required"))
        for layer in layers
            length(layer.bonds) == length(layer.gates) || throw(ArgumentError("one gate per bond required"))
            all(i -> 1 <= i < L, layer.bonds) || throw(ArgumentError("invalid open-chain bond"))
            all(i -> 1 <= i <= L, layer.measured) || throw(ArgumentError("invalid measurement site"))
            allunique(layer.measured) || throw(ArgumentError("duplicate measurement site"))
            for U in layer.gates
                size(U) == (4, 4) && all(isfinite, U) && isapprox(U' * U, I; atol=1e-11, rtol=1e-11) ||
                    throw(ArgumentError("gates must be finite 4x4 unitaries"))
            end
        end
        new(Int(L), layers)
    end
end

function haar_u2(rng)
    alpha, psi, chi = 2pi .* rand(rng, 3)
    theta = asin(sqrt(rand(rng)))
    exp(im * alpha) * [exp(im * psi)*cos(theta) exp(im * chi)*sin(theta);
                      -exp(-im * chi)*sin(theta) exp(-im * psi)*cos(theta)]
end

function dual_gate(rng)
    phase, J = 2pi .* rand(rng, 2)
    up, um, vp, vm = (haar_u2(rng) for _ in 1:4)
    # Closed form of exp[-i pi/4 (XX + YY) - i J ZZ], matching the archive.
    V = ComplexF64[exp(-im*J) 0 0 0;
                  0 0 -im*exp(im*J) 0;
                  0 -im*exp(im*J) 0 0;
                  0 0 0 exp(-im*J)]
    exp(im * phase) * kron(up, um) * V * kron(vp, vm)
end

function random_circuit(rng::AbstractRNG; L::Integer=6, p::Real=0.14,
                        layers::Integer=16, dephase::Bool=true)
    iseven(L) && 2 <= L <= 20 || throw(ArgumentError("L must be even, from 2 to 20"))
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must be in [0,1]"))
    1 <= layers <= 10000 || throw(ArgumentError("layers must be from 1 to 10000"))
    plan = Layer[]
    for t in 1:layers
        bonds = collect(isodd(t) ? (1:2:L-1) : (2:2:L-2))
        gates = [dual_gate(rng) for _ in bonds]
        sites = isodd(t) ? (1:L) : (2:L-1)
        measured = [i for i in sites if rand(rng) < p]
        push!(plan, Layer(bonds, gates, measured; dephase=dephase && iseven(t)))
    end
    Circuit(L, plan)
end

function initial_state(rng::AbstractRNG, L::Integer; kind::Symbol=:haar)
    iseven(L) && 2 <= L <= 20 || throw(ArgumentError("invalid L"))
    kind in (:haar, :product) || throw(ArgumentError("kind must be haar or product"))
    if kind == :haar
        return normalize!(randn(rng, ComplexF64, 2^L))
    end
    psi = zeros(ComplexF64, 2^L)
    psi[rand(rng, 1:length(psi))] = 1
    psi
end

function checked_state(c, initial)
    initial isa AbstractVector && length(initial) == 2^c.L && all(isfinite, initial) &&
        isapprox(norm(initial), 1; atol=1e-11, rtol=1e-11) || throw(ArgumentError("initial state must be normalized"))
    ComplexF64.(initial)
end

function check_record(c, record)
    length(record) == length(c.layers) || throw(ArgumentError("one outcome list per layer required"))
    for (layer, outcomes) in zip(c.layers, record)
        length(layer.measured) == length(outcomes) && all(x -> x in (0, 1), outcomes) ||
            throw(ArgumentError("record must contain one bit per measured site"))
    end
    nothing
end

function apply_gate!(psi, U, site)
    stride = 1 << (site - 1)
    @inbounds for base in 0:4stride:length(psi)-1, offset in 1:stride
        a = base + offset
        b, c, d = a + stride, a + 2stride, a + 3stride
        v1, v2, v3, v4 = psi[a], psi[b], psi[c], psi[d]
        psi[a] = U[1,1]*v1 + U[1,2]*v2 + U[1,3]*v3 + U[1,4]*v4
        psi[b] = U[2,1]*v1 + U[2,2]*v2 + U[2,3]*v3 + U[2,4]*v4
        psi[c] = U[3,1]*v1 + U[3,2]*v2 + U[3,3]*v3 + U[3,4]*v4
        psi[d] = U[4,1]*v1 + U[4,2]*v2 + U[4,3]*v3 + U[4,4]*v4
    end
    psi
end

function project!(psi, site, outcome)
    mask = 1 << (site - 1)
    @inbounds for k in eachindex(psi)
        ((k - 1) & mask != 0) == Bool(outcome) || (psi[k] = 0)
    end
    q = sum(abs2, psi)
    q > 0 && (psi ./= sqrt(q))
    clamp(q, 0.0, 1.0)
end

function phase_flip!(psi)
    @views psi[length(psi)÷2+1:end] .*= -1
    psi
end

"""Draw the observed record with its Born marginal; hidden boundary jumps are discarded."""
function sample_record(rng::AbstractRNG, c::Circuit, initial)
    psi = checked_state(c, initial)
    record = Vector{Int}[]
    for layer in c.layers
        for (i, U) in zip(layer.bonds, layer.gates)
            apply_gate!(psi, U, i)
        end
        outcomes = Int[]
        for site in layer.measured
            mask = 1 << (site - 1)
            q0 = sum(abs2(psi[k]) for k in eachindex(psi) if (k-1) & mask == 0)
            outcome = Int(rand(rng) >= q0)
            project!(psi, site, outcome) > 0 || error("sampled an impossible outcome")
            push!(outcomes, outcome)
        end
        push!(record, outcomes)
        layer.dephase && rand(rng, Bool) && phase_flip!(psi)
    end
    record
end

function path_logprob_unchecked(c, initial, record, jumps)
    psi = copy(initial)
    values = zeros(length(c.layers))
    logp, jump = 0.0, 0
    for (t, layer) in enumerate(c.layers)
        for (i, U) in zip(layer.bonds, layer.gates)
            apply_gate!(psi, U, i)
        end
        for (site, outcome) in zip(layer.measured, record[t])
            q = project!(psi, site, outcome)
            if q == 0
                values[t:end] .= -Inf
                return values
            end
            logp += log(q)
        end
        values[t] = logp
        if layer.dephase
            jump += 1
            jumps[jump] && phase_flip!(psi)
        end
    end
    values
end

"""Log P(record prefix | hidden I/Z history), excluding the uniform history weight."""
function path_logprob(c::Circuit, initial, record, jumps)
    psi = checked_state(c, initial)
    check_record(c, record)
    length(jumps) == count(l -> l.dephase, c.layers) && all(x -> x in (0, 1), jumps) ||
        throw(ArgumentError("one binary jump per dephasing event required"))
    path_logprob_unchecked(c, psi, record, Bool.(jumps))
end

function embed_gate(U, site, L)
    left = spdiagm(0 => ones(ComplexF64, 2^(L-site-1)))
    right = spdiagm(0 => ones(ComplexF64, 2^(site-1)))
    kron(left, sparse(U), right)
end

"""Full 2^L by 2^L conditional density matrix, with exact boundary dephasing."""
function density_record(c::Circuit, initial, record)
    c.L <= 10 || throw(ArgumentError("full density-matrix reference is limited to L <= 10"))
    psi = checked_state(c, initial)
    check_record(c, record)
    rho = psi * psi'
    values = zeros(length(c.layers))
    logp = 0.0
    for (t, layer) in enumerate(c.layers)
        for (site, U) in zip(layer.bonds, layer.gates)
            A = embed_gate(U, site, c.L)
            rho = A * rho * A'
        end
        for (site, outcome) in zip(layer.measured, record[t])
            mask = 1 << (site - 1)
            rejected = [k for k in axes(rho, 1) if (((k-1) & mask) != 0) != Bool(outcome)]
            rho[rejected, :] .= 0
            rho[:, rejected] .= 0
            q = real(tr(rho))
            if q <= 0
                values[t:end] .= -Inf
                return (; logprob=values, rho)
            end
            rho ./= q
            logp += log(clamp(q, 0.0, 1.0))
        end
        if layer.dephase
            half = size(rho, 1) ÷ 2
            rho[1:half, half+1:end] .= 0
            rho[half+1:end, 1:half] .= 0
        end
        values[t] = logp
    end
    (; logprob=values, rho)
end

function logadd(a, b)
    a == -Inf && return b
    b == -Inf && return a
    max(a, b) + log1p(exp(-abs(a-b)))
end

"""Exact uniform sum over every hidden I/Z history, for small circuits."""
function enumerate_histories(c::Circuit, initial, record)
    psi = checked_state(c, initial)
    check_record(c, record)
    n = count(l -> l.dephase, c.layers)
    n <= 16 || throw(ArgumentError("exhaustive enumeration is limited to 16 dephasing events"))
    total = fill(-Inf, length(c.layers))
    for index in 0:2^n-1
        jumps = [Bool((index >> (i-1)) & 1) for i in 1:n]
        values = path_logprob_unchecked(c, psi, record, jumps)
        total .= logadd.(total, values)
    end
    total .- n * log(2)
end

"""Streaming Monte Carlo estimate of record probabilities at cumulative sample counts.

Columns correspond to counts; rows to physical gate/measurement layers.
relative_se estimates the standard error of the probability divided by its mean.
ess reports importance-weight effective sample size. Zero estimates remain -Inf.
"""
function monte_carlo(rng::AbstractRNG, c::Circuit, initial, record;
                     counts=[100, 300, 400, 500, 600])
    psi = checked_state(c, initial)
    check_record(c, record)
    !isempty(counts) && all(x -> x isa Integer && x >= 2, counts) &&
        issorted(counts) && allunique(counts) || throw(ArgumentError("counts must be increasing integers >= 2"))
    ns = Int.(counts)
    T = length(c.layers)
    njumps = count(l -> l.dephase, c.layers)
    sums, squares = fill(-Inf, T), fill(-Inf, T)
    logprob, relative_se, ess = (zeros(T, length(ns)) for _ in 1:3)
    checkpoint = 1
    for n in 1:last(ns)
        values = path_logprob_unchecked(c, psi, record, rand(rng, Bool, njumps))
        sums .= logadd.(sums, values)
        squares .= logadd.(squares, 2 .* values)
        if n == ns[checkpoint]
            for t in 1:T
                logprob[t, checkpoint] = sums[t] - log(n)
                if sums[t] == -Inf
                    relative_se[t, checkpoint], ess[t, checkpoint] = Inf, 0.0
                else
                    ratio = log(n) + squares[t] - 2sums[t]
                    relative_se[t, checkpoint] = sqrt(max(expm1(ratio), 0.0) / (n-1))
                    ess[t, checkpoint] = min(n, exp(2sums[t] - squares[t]))
                end
            end
            checkpoint += 1
        end
    end
    (; counts=ns, logprob, relative_se, ess)
end

function free_energy_density(logprob::AbstractVector, L::Integer; burnin::Integer=0)
    L > 0 && 0 <= burnin < length(logprob) || throw(ArgumentError("invalid L or burnin"))
    all(isfinite, logprob) || throw(ArgumentError("non-finite record likelihood; increase sampling or inspect the record"))
    -(last(logprob) - (burnin == 0 ? 0.0 : logprob[burnin])) / (L * (length(logprob)-burnin))
end

"""Fit ensemble free energies to a polynomial in 1/tau and return its tau=infinity intercept."""
function extrapolate_samples(counts, values; degree::Integer=4, min_samples::Integer=300)
    length(counts) == length(values) && all(isfinite, values) || throw(ArgumentError("invalid fit data"))
    all(x -> x isa Integer && x > 0, counts) && allunique(counts) || throw(ArgumentError("invalid sample counts"))
    degree in 1:4 && min_samples > 0 || throw(ArgumentError("degree must be 1:4 and min_samples positive"))
    selected = findall(>=(min_samples), counts)
    length(selected) >= degree + 2 || throw(ArgumentError("too few sample counts for requested fit"))
    scale = minimum(counts[selected])
    x = scale ./ counts[selected]
    X = hcat([x.^k for k in 0:degree]...)
    coefficients = qr(X) \ values[selected]
    (; intercept=coefficients[1], coefficients, scale, condition=cond(X),
       residual_rms=sqrt(mean(abs2, X * coefficients - values[selected])))
end

"""Resample complete circuit/record curves to retain correlations between tau checkpoints."""
function bootstrap_extrapolation(rng::AbstractRNG, counts, curves::AbstractMatrix;
                                 degree=4, min_samples=300, replicates::Integer=500)
    size(curves, 1) >= 2 && size(curves, 2) == length(counts) && all(isfinite, curves) ||
        throw(ArgumentError("curves must have >= 2 independent rows and one column per count"))
    replicates >= 20 || throw(ArgumentError("at least 20 bootstrap replicates required"))
    estimate = extrapolate_samples(counts, vec(mean(curves; dims=1)); degree, min_samples)
    intercepts = zeros(replicates)
    for b in 1:replicates
        ids = rand(rng, 1:size(curves, 1), size(curves, 1))
        intercepts[b] = extrapolate_samples(counts, vec(mean(curves[ids, :]; dims=1)); degree, min_samples).intercept
    end
    (; estimate..., stderr=std(intercepts), interval95=quantile(intercepts, [0.025, 0.975]))
end

"""Fit L*(f_dissipative - f_open) = surface_difference + pi*h/L."""
function fit_boundary_dimension(sizes, deltaf; Lmin=minimum(sizes))
    length(sizes) == length(deltaf) && all(isfinite, deltaf) &&
        all(x -> x isa Integer && x > 0, sizes) && allunique(sizes) || throw(ArgumentError("invalid finite-size data"))
    ids = findall(>=(Lmin), sizes)
    length(ids) >= 3 || throw(ArgumentError("at least three sizes required"))
    L = Float64.(sizes[ids])
    X, y = hcat(ones(length(L)), 1 ./ L), L .* deltaf[ids]
    coef = qr(X) \ y
    (; surface_difference=coef[1], h=coef[2]/pi, residual_rms=sqrt(mean(abs2, X * coef-y)))
end

end
