module TestExamplesThreaded

using Test
using Trixi

const TEST_DIR = joinpath((dirname ∘ dirname ∘ pathof)(Trixi), "test")

include(joinpath(TEST_DIR, "test_trixi.jl"))

EXAMPLES_DIR = examples_dir()

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
Trixi.mpi_isroot() && isdir(outdir) && rm(outdir, recursive = true)
Trixi.MPI.Barrier(Trixi.mpi_comm())

@testset "Threaded tests" begin
#! format: noindent

@testset "TreeMesh" begin
    @trixi_testset "elixir_advection_restart.jl" begin
        elixir = joinpath(EXAMPLES_DIR, "tree_2d_dgsem",
                          "elixir_advection_extended.jl")
        mpi_isroot() && println("═"^100)
        mpi_isroot() && println(elixir)
        trixi_include(@__MODULE__, elixir, tspan = (0.0, 10.0))
        l2_expected, linf_expected = analysis_callback(sol)

        elixir = joinpath(EXAMPLES_DIR, "tree_2d_dgsem",
                          "elixir_advection_restart.jl")
        mpi_isroot() && println("═"^100)
        mpi_isroot() && println(elixir)
   end
end
end
end # module