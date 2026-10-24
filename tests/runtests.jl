using Test, LinearAlgebra, Random
include("../examples/quickstart.jl")

@testset "Original open-boundary kernels" begin
    for seed in 1:10
        Random.seed!(seed)
        U = OpenBoundary.rangate()
        @test U' * U ≈ I atol=1e-12
        reshuffled = reshape(permutedims(reshape(U, 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
        @test reshuffled' * reshuffled ≈ I atol=1e-12
    end
    values = OpenBoundary.trajectory(; L=4, steps=24, seed=3)
    @test values == OpenBoundary.trajectory(; L=4, steps=24, seed=3)
    @test all(isfinite, values)
    @test minimum(diff(values)) >= -1e-12
    @test OpenBoundary.trajectory(; p=0, steps=12) ≈ zeros(12) atol=1e-12
    OpenBoundary.trajectory(; L=4, steps=1)
    psi = normalize(randn(ComplexF64, 16))
    for i in (1, 2, 3)
        Random.seed!(22)
        U = OpenBoundary.rangate()
        dense = kron(Matrix{ComplexF64}(I, 2^(4-i-1), 2^(4-i-1)), U,
                     Matrix{ComplexF64}(I, 2^(i-1), 2^(i-1)))
        Random.seed!(22)
        @test OpenBoundary.gate(copy(psi), i) ≈ dense * psi atol=1e-12
    end
    @test OpenBoundary.gate(copy(psi), 4) == psi
end

@testset "Measurement record and density-matrix reference" begin
    OpenBoundary.trajectory(; L=2, p=1, steps=1)
    global OpenBoundary.t = 1
    global OpenBoundary.FET = zeros(1)
    psi = ComplexF64[1, 0, 0, 1] / sqrt(2)
    after = OpenBoundary.even_projmsm(psi)
    @test norm(after) ≈ 1
    @test OpenBoundary.FET[1] ≈ log(2)
    @test count(x -> abs(x) > 1e-12, after) == 1
    # Both projectors give probability 1/2 for the Bell state.
    rho = ComplexF64[1, 0, 0, 1] * ComplexF64[1, 0, 0, 1]' / 2
    P0 = Diagonal([1, 0, 1, 0])
    P1 = I - P0
    @test real(tr(P0 * rho * P0)) ≈ 0.5
    @test real(tr(P1 * rho * P1)) ≈ 0.5
end

@testset "Input validation and safe outputs" begin
    for kwargs in ((; L=3), (; p=NaN), (; p=-0.1), (; steps=0), (; seed=-1))
        @test_throws ArgumentError OpenBoundary.trajectory(; kwargs...)
    end
    mktempdir() do dir
        out = joinpath(dir, "demo")
        main([out])
        @test isfile(joinpath(out, "trace.csv"))
        saved = read(joinpath(out, "trace.csv"))
        @test_throws ArgumentError main([out])
        @test saved == read(joinpath(out, "trace.csv"))
    end
end
