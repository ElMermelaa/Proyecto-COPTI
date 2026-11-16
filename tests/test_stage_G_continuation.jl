using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "scripts", "run_stage_G_presentation_videos.jl"))

@testset "Etapa G — continuación física después de T" begin
    for (name, kind) in (("etapa_E_sync.jld2", "synchronization"),
                         ("etapa_E_desync.jld2", "desynchronization"))
        file = joinpath(DATA_G, name)
        @test isfile(file)
        old = read_stage_e(file, kind)
        continuation = extend_stage_e(file, kind;
                         post_time=4.0, post_control=:hold)
        @test continuation.n > old.n
        @test continuation.times[old.n] == old.times[end]
        @test isapprox(last(continuation.times), last(old.times)+4.0; atol=1e-9)
        @test continuation.X[:,1:old.n] == old.X
        @test all(diff(continuation.times) .> 0)
        # Comprueba que DESPUÉS de T existen estados nuevos, NO cuadros repetidos.
        @test norm(continuation.X[:,end] - old.X[:,end]) > 1e-6
        @test isapprox(continuation.E[end],
                       synchronization_error(continuation.X[:,end]); atol=1e-12)
        if kind == "synchronization"
            @test length(continuation.control) == continuation.n
        else
            @test size(continuation.control) == (old.N, continuation.n)
        end
    end
end
println("Etapa G: se verificó la continuación real en ambas simulaciones.")
