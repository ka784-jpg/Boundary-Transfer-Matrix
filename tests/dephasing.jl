include("../examples/dephasing.jl")
const BD = DephasingExample.BoundaryDephasing
BLAS.set_num_threads(1)

@testset "Dual-unitary gates and statevector ordering" begin
    rng = Xoshiro(711)
    for L in (2, 4, 6), i in 1:L-1
        U = BD.dual_gate(rng)
        psi = BD.initial_state(rng, L)
        @test U' * U ≈ I atol=1e-12
        reshuffled = reshape(permutedims(reshape(U, 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
        @test reshuffled' * reshuffled ≈ I atol=1e-12
        dense = kron(Matrix{ComplexF64}(I, 2^(L-i-1), 2^(L-i-1)), U,
                     Matrix{ComplexF64}(I, 2^(i-1), 2^(i-1)))
        @test BD.apply_gate!(copy(psi), U, i) ≈ dense * psi atol=1e-12
    end
end

@testset "Analytical channel, record probabilities, and hidden entropy" begin
    H = ComplexF64[1 1; 1 -1] / sqrt(2)
    U = kron(H, Matrix{ComplexF64}(I, 2, 2))
    psi = ComplexF64[1, 0, 0, 0]
    c = BD.Circuit(2, [BD.Layer([1], [U], Int[]; dephase=true), BD.Layer([1], [U], [2])])
    probabilities = [exp(BD.density_record(c, psi, [Int[], [bit]]).logprob[end]) for bit in 0:1]
    @test probabilities ≈ [0.5, 0.5] atol=1e-13
    @test sum(probabilities) ≈ 1
    @test -sum(probabilities .* log.(probabilities)) ≈ log(2)
    @test BD.enumerate_histories(c, psi, [Int[], [0]]) ≈ [0, -log(2)] atol=1e-13
    mc = BD.monte_carlo(Xoshiro(22), c, psi, [Int[], [0]]; counts=[64, 4096])
    @test abs(exp(mc.logprob[end, end]) - 0.5) < 0.03
    @test isfinite(mc.logprob[end, end])
    zeros_count = count(_ -> last(BD.sample_record(Random.default_rng(), c, psi)) == [0], 1:4096)
    @test abs(zeros_count/4096 - 0.5) < 0.04
    no_observation = BD.Circuit(2, [BD.Layer([1], [U], Int[]; dephase=true)])
    rho = BD.density_record(no_observation, psi, [Int[]]).rho
    pure = U * psi
    Z = Diagonal([1, 1, -1, -1])
    @test rho ≈ (pure * pure' + Z * pure * pure' * Z)/2 atol=1e-13
    @test BD.monte_carlo(Xoshiro(3), no_observation, psi, [Int[]]; counts=[16]).logprob[1] ≈ 0 atol=1e-13
end

@testset "Full density matrix versus every hidden history" begin
    for L in (2, 4, 6, 8), kind in (:haar, :product), p in (0.0, 0.14, 0.7, 1.0), dissipative in (false, true)
        rng = Xoshiro(100L + round(Int, 100p) + Int(dissipative))
        c = BD.random_circuit(rng; L, p, layers=6, dephase=dissipative)
        psi = BD.initial_state(rng, L; kind)
        record = BD.sample_record(rng, c, psi)
        exact = BD.density_record(c, psi, record)
        @test BD.enumerate_histories(c, psi, record) ≈ exact.logprob atol=2e-11
        @test tr(exact.rho) ≈ 1 atol=1e-12
        @test exact.rho ≈ exact.rho' atol=1e-12
        @test all(diff(exact.logprob) .<= 1e-12)
    end
end

@testset "Monte Carlo convergence, diagnostics, and repeatability" begin
    for seed in 1:8
        rng = Xoshiro(seed)
        c = BD.random_circuit(rng; L=4, p=0.5, layers=8)
        psi = BD.initial_state(rng, 4)
        record = BD.sample_record(rng, c, psi)
        exact = exp.(BD.density_record(c, psi, record).logprob)
        mc = BD.monte_carlo(Xoshiro(123seed), c, psi, record; counts=[64, 512, 4096])
        estimate = exp.(mc.logprob[:, end])
        @test all(abs.(estimate - exact) .<= 8 .* estimate .* mc.relative_se[:, end] .+ 1e-12)
        @test all(1 - 1e-10 .<= mc.ess .<= reshape(mc.counts, 1, :) .+ 1e-10)
        @test all(mc.logprob .<= 1e-12)
        @test mc.logprob == BD.monte_carlo(Xoshiro(123seed), c, psi, record; counts=mc.counts).logprob
    end
end

@testset "Log-domain stability and impossible records" begin
    H = ComplexF64[1 1; 1 -1] / sqrt(2)
    U = kron(H, Matrix{ComplexF64}(I, 2, 2))
    psi = ComplexF64[1, 0, 0, 0]
    c = BD.Circuit(2, [BD.Layer([1], [U], [2]; dephase=true) for _ in 1:1200])
    record = [[0] for _ in c.layers]
    exact = BD.density_record(c, psi, record)
    mc = BD.monte_carlo(Xoshiro(4), c, psi, record; counts=[2])
    @test exact.logprob[end] ≈ -1200log(2) atol=1e-9
    @test mc.logprob[end] ≈ exact.logprob[end] atol=1e-9
    @test exp(mc.logprob[end]) == 0
    impossible = BD.Circuit(2, [BD.Layer(Int[], Matrix{ComplexF64}[], [1])])
    @test BD.density_record(impossible, psi, [[1]]).logprob[end] == -Inf
    @test BD.monte_carlo(Xoshiro(1), impossible, psi, [[1]]; counts=[2]).logprob[end] == -Inf
    @test_throws ArgumentError BD.free_energy_density([-Inf], 2)
end

@testset "Extrapolation, bootstrap, and finite-size fitting" begin
    counts = collect(300:50:600)
    x = 300 ./ counts
    values = @. 0.075 + 0.003x - 0.0005x^2 + 0.0002x^3 - 0.00002x^4
    fit = BD.extrapolate_samples(counts, values)
    @test fit.intercept ≈ 0.075 atol=1e-11
    @test fit.residual_rms < 1e-12
    curves = repeat(values', 4, 1) .+ [-0.003, -0.001, 0.001, 0.003]
    bootstrap = BD.bootstrap_extrapolation(Xoshiro(3), counts, curves; replicates=100)
    @test bootstrap.intercept ≈ 0.075 atol=1e-11
    @test bootstrap.stderr > 0
    @test bootstrap.interval95[1] < 0.075 < bootstrap.interval95[2]
    sizes = [8, 10, 12, 14, 16]
    delta = @. 0.004 / sizes + pi * 0.05 / sizes^2
    @test BD.fit_boundary_dimension(sizes, delta).h ≈ 0.05 atol=1e-12
    @test BD.fit_boundary_dimension(sizes, delta; Lmin=12).surface_difference ≈ 0.004 atol=1e-12
    @test BD.free_energy_density([0.0, -1.0, -3.0], 2; burnin=1) == 0.75
end

@testset "Dephasing input contracts" begin
    rng = Xoshiro(2)
    c = BD.random_circuit(rng; L=2, layers=2)
    psi = BD.initial_state(rng, 2)
    record = BD.sample_record(rng, c, psi)
    @test_throws ArgumentError BD.random_circuit(rng; L=3)
    @test_throws ArgumentError BD.random_circuit(rng; p=NaN)
    @test_throws ArgumentError BD.initial_state(rng, 2; kind=:unknown)
    @test_throws ArgumentError BD.density_record(c, 2psi, record)
    @test_throws ArgumentError BD.density_record(c, psi, [[2]])
    @test_throws ArgumentError BD.path_logprob(c, psi, record, Bool[])
    for counts in ([1], [4, 2], [2, 2], Int[])
        @test_throws ArgumentError BD.monte_carlo(rng, c, psi, record; counts)
    end
    @test_throws ArgumentError BD.extrapolate_samples([300, 400], [1.0, 2.0])
    @test_throws ArgumentError BD.fit_boundary_dimension([4, 6], [0.1, 0.2])
    @test_throws ArgumentError BD.free_energy_density([0.0], 2; burnin=1)
    mktempdir() do dir
        @test_throws ArgumentError DephasingExample.run(joinpath(@__DIR__, "..", "examples", "dephasing.toml"), dir)
    end
end
