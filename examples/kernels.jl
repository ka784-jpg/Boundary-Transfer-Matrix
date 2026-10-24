module OpenBoundary
using LinearAlgebra, Random, Statistics, Distributions, DelimitedFiles, Dates

# Supply Normal explicitly; the archive relied on its surrounding environment.
function definitions_only(ex)
    ex isa Expr || return ex
    ex.head in (:using, :import) && return nothing
    ex.head == :call && ex.args[1] == :runy && return nothing
    return ex
end
Base.include(definitions_only, @__MODULE__, joinpath(@__DIR__, "..", "research", "open_boundary_haar.jl"))

function trajectory(; L=4, p=0.14, steps=48, seed=2026)
    L isa Integer && iseven(L) && 2 <= L <= 10 || throw(ArgumentError("L must be even, 2 through 10"))
    p isa Real && isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must be in [0,1]"))
    steps isa Integer && 1 <= steps <= 320 || throw(ArgumentError("steps must be 1 through 320"))
    seed isa Integer && 0 <= seed <= typemax(Int32) || throw(ArgumentError("invalid seed"))
    global nspin = L
    global prob = Float64(p)
    global FET = zeros(steps)
    Random.seed!(seed)
    psi = randn(ComplexF64, 2^L)
    psi /= norm(psi)
    for layer in 1:steps
        global t = layer
        bonds = isodd(layer) ? (1:2:L-1) : (2:2:L-2)
        for i in bonds
            psi = gate(psi, i)
        end
        psi = isodd(layer) ? even_projmsm(psi) : odd_projmsm(psi)
        isfinite(norm(psi)) && isapprox(norm(psi), 1; atol=1e-10) || error("non-normalized state")
    end
    all(isfinite, FET) || error("non-finite measurement-record surprisal")
    return copy(FET)
end
end
