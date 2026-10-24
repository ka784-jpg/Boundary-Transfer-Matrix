module DephasingExample
using Random, Statistics, LinearAlgebra, DelimitedFiles, TOML, SHA
include("../src/BoundaryDephasing.jl")
using .BoundaryDephasing

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join(row, ','))
        end
    end
end

function run(config_path, output_path)
    config = TOML.parsefile(config_path)
    allowed = Set(["sizes", "p", "layers", "records_per_state", "initial_states", "counts",
                   "seed", "burnin", "density_reference", "exhaustive", "bootstrap_replicates"])
    isempty(setdiff(keys(config), allowed)) || throw(ArgumentError("unrecognized configuration key"))
    sizes = config["sizes"]
    allunique(sizes) && !isempty(sizes) || throw(ArgumentError("sizes must be nonempty and distinct"))
    nrecords = config["records_per_state"]
    nrecords isa Integer && nrecords >= 2 || throw(ArgumentError("at least two records per initial state required"))
    kinds = Symbol.(config["initial_states"])
    !isempty(kinds) && allunique(kinds) && all(k -> k in (:haar, :product), kinds) ||
        throw(ArgumentError("invalid initial_states"))
    seed = config["seed"]
    seed isa Integer && 0 <= seed <= typemax(Int32) || throw(ArgumentError("invalid seed"))
    layers, burnin, counts = config["layers"], config["burnin"], config["counts"]
    burnin isa Integer && 0 <= burnin < layers || throw(ArgumentError("invalid burnin"))
    reference, exhaustive = config["density_reference"], config["exhaustive"]
    reference isa Bool && exhaustive isa Bool || throw(ArgumentError("reference flags must be boolean"))
    reference && any(>(10), sizes) && throw(ArgumentError("density reference requires L <= 10"))
    exhaustive && layers > 32 && throw(ArgumentError("exhaustive check requires <= 32 layers"))
    output = abspath(output_path)
    ispath(output) && throw(ArgumentError("output already exists: $output"))
    mkpath(dirname(output))
    stage = mktempdir(dirname(output))
    try
        trajectories, curves, fits = Any[], Any[], Any[]
        max_enumeration_error = 0.0
        for L in sizes, (kind_index, kind) in enumerate(kinds)
            samples = zeros(nrecords, length(counts))
            exact = reference ? zeros(nrecords) : nothing
            for record_id in 1:nrecords
                case_seed = seed + 100000L + 10000kind_index + record_id
                rng = Xoshiro(case_seed)
                circuit = random_circuit(rng; L, p=config["p"], layers)
                psi = initial_state(rng, L; kind)
                record = sample_record(rng, circuit, psi)
                # A separate stream avoids correlating likelihood estimation with record generation.
                mc = monte_carlo(Xoshiro(case_seed + 1000000000), circuit, psi, record; counts)
                dm = reference ? density_record(circuit, psi, record).logprob : nothing
                enumerated = exhaustive ? enumerate_histories(circuit, psi, record) : nothing
                if reference && exhaustive
                    discrepancy = maximum(abs, enumerated - dm)
                    max_enumeration_error = max(max_enumeration_error, discrepancy)
                    discrepancy < 5e-10 || error("exhaustive statevector / density-matrix mismatch")
                end
                reference && (exact[record_id] = free_energy_density(dm, L; burnin))
                for (j, tau) in enumerate(counts)
                    samples[record_id, j] = free_energy_density(mc.logprob[:, j], L; burnin)
                    push!(trajectories, (L, kind, record_id, case_seed, tau, samples[record_id, j],
                        reference ? exact[record_id] : "", mc.logprob[end, j],
                        reference ? dm[end] : "", mc.relative_se[end, j], mc.ess[end, j],
                        reference && exhaustive ? maximum(abs, enumerated - dm) : ""))
                end
            end
            for (j, tau) in enumerate(counts)
                paired = reference ? samples[:, j] - exact : nothing
                push!(curves, (L, kind, tau, mean(samples[:, j]), std(samples[:, j])/sqrt(nrecords),
                    reference ? mean(exact) : "", reference ? std(exact)/sqrt(nrecords) : "",
                    reference ? mean(paired) : "", reference ? std(paired)/sqrt(nrecords) : ""))
            end
            if count(>=(300), counts) >= 6
                fit = bootstrap_extrapolation(Xoshiro(seed + 100L + kind_index), counts, samples;
                                              replicates=config["bootstrap_replicates"])
                push!(fits, (L, kind, fit.intercept, fit.stderr, fit.interval95...,
                             fit.condition, fit.residual_rms, reference ? mean(exact) : ""))
            end
            println("Completed L=$L, $kind: $nrecords records, $(last(counts)) hidden histories each")
        end
        write_csv(joinpath(stage, "records.csv"),
            ["L", "initial_state", "record", "seed", "tau", "mc_f", "density_f",
             "mc_logprob", "density_logprob", "relative_se", "ess", "enumeration_max_error"], trajectories)
        write_csv(joinpath(stage, "convergence.csv"),
            ["L", "initial_state", "tau", "mc_f", "mc_se", "density_f", "density_se", "paired_error", "paired_se"], curves)
        write_csv(joinpath(stage, "extrapolation.csv"),
            ["L", "initial_state", "intercept", "bootstrap_se", "ci_low", "ci_high", "condition", "residual_rms", "density_f"], fits)
        metadata = Dict("config" => config, "julia_version" => string(VERSION),
            "units" => "nats per site per gate-plus-measurement layer",
            "record_seed_rule" => "seed + 100000*L + 10000*initial_state_index + record_id",
            "inner_seed_rule" => "record_seed + 1000000000",
            "max_enumeration_logprob_error" => max_enumeration_error,
            "source_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "..", "src", "BoundaryDephasing.jl")))),
            "driver_sha256" => bytes2hex(sha256(read(@__FILE__))))
        open(joinpath(stage, "metadata.toml"), "w") do io
            TOML.print(io, metadata)
        end
        mv(stage, output)
        println("Maximum exhaustive / density log-probability difference: $max_enumeration_error")
    finally
        isdir(stage) && rm(stage; recursive=true)
    end
    output
end

function main(args)
    length(args) == 2 || throw(ArgumentError("usage: julia --project=. examples/dephasing.jl CONFIG.toml OUTPUT_DIRECTORY"))
    BLAS.set_num_threads(1)
    run(args...)
end
end

abspath(PROGRAM_FILE) == (@__FILE__) && DephasingExample.main(ARGS)
