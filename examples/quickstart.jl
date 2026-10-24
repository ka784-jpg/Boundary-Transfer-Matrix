using DelimitedFiles, TOML, SHA
include("kernels.jl")

function main(args)
    length(args) == 1 || throw(ArgumentError("usage: julia --project=. examples/quickstart.jl OUTPUT_DIRECTORY"))
    output = abspath(only(args))
    ispath(output) && throw(ArgumentError("output already exists: $output"))
    parent = dirname(output)
    mkpath(parent)
    stage = mktempdir(parent)
    try
        rows = Matrix{Float64}(undef, 0, 5)
        for L in (4, 6)
            steps = 16L
            for trajectory_id in 1:8
                seed = 2026 + 100L + trajectory_id
                values = OpenBoundary.trajectory(; L, p=0.14, steps, seed)
                rows = vcat(rows, hcat(fill(L, steps), fill(trajectory_id, steps),
                            1:steps, values, values ./ (L .* (1:steps))))
            end
        end
        open(joinpath(stage, "trace.csv"), "w") do io
            println(io, "L,trajectory,layer,surprisal_nats,surprisal_per_site_layer")
            writedlm(io, rows, ',')
        end
        open(joinpath(stage, "metadata.toml"), "w") do io
            TOML.print(io, Dict("p" => 0.14, "trajectories_per_size" => 8,
                "sizes" => [4, 6], "layers_per_site" => 16,
                "seed_rule" => "2026 + 100*L + trajectory_id", "julia_version" => string(VERSION),
                "description" => "Small open-boundary example, no dephasing channel or CFT fit",
                "source_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "..", "research", "open_boundary_haar.jl"))))))
        end
        mv(stage, output)
    finally
        isdir(stage) && rm(stage; recursive=true)
    end
    println("Saved open-boundary example to $output")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
