# In this file we setup the testing framework for our package
# The Julia documentation for testing can be found at:
# https://docs.julialang.org/en/v1/stdlib/Test/#Base.runtests
using Test
using MyPkg

# Define the tests

@testset "MyPkg Tests" begin
    # Example test
    @test foo(2) == 4
    @test foo2(2) == 8
    # Add more tests here
end

# Run some tests from different files if you defined them. Especially useful for large packages

function runtests()
    testdir = pwd()
    istest(f) = endswith(f, ".jl") && startswith(basename(f), "test_")
    testfiles = sort(
        filter(
            istest,
            vcat([joinpath.(root, files) for (root, dirs, files) in walkdir(testdir)]...),
        ),
    )
    nfail = 0
    printstyled("Testing package MyPkg.jl\n"; bold=true, color=:white)

    for f in testfiles
        println("")
        println("Running tests from $f")

        try
            run(`$(Base.julia_cmd()) -O3 --startup-file=no --check-bounds=no $(joinpath(testdir, f))`)
        catch ex
            nfail += 1
        end
    end

    return nfail
end

runtests()
