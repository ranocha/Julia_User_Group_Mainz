# # Essential Julia packages and debugging

# First, you need to install Julia. I use juliaup:
#
# - https://julialang.org/downloads/
#
# To develop some Julia code, you need an editor and a way to run code.
# I prefer Visual Studio Code with the Julia extension as editor since
# it comes with a nice integrated development environment (IDE).
# You can also run Julia via VS code. I will likely do so for this session.
# However, I prefer having a separate terminal open when developing code.
#
# - VS code: https://code.visualstudio.com/
# - Julia extension: https://www.julia-vscode.org/
#
# Here is a summary of the Julia workflow we will cover:
#
# 0. Install Revise.jl, e.g., in you default Julia environment
#    ```julia
#    using Pkg
#    Pkg.activate() # activate the default environment/project
#    Pkg.add("Revise")
#    ```
# 1. Create a new folder for your project
# 2. Open the Julia REPL in this folder
# 3. Create a Julia environment and add packages you want to use
#    ```julia
#    using Pkg
#    Pkg.activate(".") # use the current folder as project
#    Pkg.add(["Plots", "LaTeXStrings"])
#    ```
# 4. Create a new Julia file for your project and let Revise track it
#    ```julia
#    using Revise
#    includet("julia_workflow.jl")
#    ```
# 5. Write code in your Julia file, execute it in the REPL, update, and
#    iterate until convergence.
# 6. Write up instructions how to reproduce your results in a `README.md`
#    file (including the Julia version you used)
# 7. Distribute (an archive of) the project folder including your Julia
#    source files, the TOML files (`Project.toml`, `Manifest.toml`) and
#    mention the Julia version you used.



# Activate the Julia project
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using LinearAlgebra
using SparseArrays
using Plots, LaTeXStrings


# Write some code
function some_linear_algebra(n)
    A = randn(n, n)
    b = randn(n)
    x = A \ b
    @info "Computation done" A b x norm(A * x - b)
    return x
end



"""
    some_sparse_stuff(n)

This is a docstring for the function `some_sparse_stuff`. This function is
written as a simple demonstration how to set up a sparse matrix of size
`n` × `n` in Julia.
"""
function some_sparse_stuff(n)
    A = spdiagm(0 => ones(n), 1 => -ones(n - 1))
end



"""
    some_plotting(n)

This is another docstring. The function `some_plotting` demonstrates how to
create a simple plot of `n` data points using the Julia package Plots.jl.
"""
function some_plotting(n)
    x = range(0.0, 1.0, length = n)
    u = sin.(π .* x) # broadcasting, can also be written as u = @. sin(π * x)
    plot(x, u, xguide = L"x", yguide = L"u", label = "sin")
end



"""
    some_conditionals_and_loops(n)

The function `some_conditionals_and_loops` demonstrates how to write
loops and conditional expressions in Julia. It accepts an integer `n`
and performs some operations, printing stuff to the REPL - and returns
`nothing`.
"""
function some_conditionals_and_loops(n::Integer)
    if n < 10
        println("n = $n < 10")
        1 + 1 # nothing is printed automatically - no need to append ;
        println(2)
    end
    println()

    for i = 1:n
        println("i = ", i)
    end
    println()

    while n^2 > 1
        println("Current value of n: ", n)
        n = n ÷ 2
    end

    return nothing
end



# ## Revise.jl
#
# - https://github.com/timholy/Revise.jl
# - import Pkg; Pkg.add("Revise")
# - using Revise; includet("path/to/your/file.jl")

function demonstrate_revise(n::Integer)
    while n^2 > 1
        println("Now, n = ", n)
        # @show n
        # @info "Test" n @__LINE__
        n = n ÷ 2
    end

    return n
end



# ## Debugging
#
# I often just use plain print statements for debugging - enhanced with
# `@show` and `@info` etc. For more complicated problems, you can either
# use a debugger or the following trick. Debuggers include
#
# - https://github.com/JuliaDebug/Infiltrator.jl
# - https://github.com/JuliaDebug/Debugger.jl (included in VS Code)
#
# Here, we demonstrate a simple, manual approach.

using OrdinaryDiffEq

function demonstrate_debugging()
    rhs(u, p, t) = -u
    u0 = [1.0, 2.0]
    tspan = (0.0, 1.0)
    ode = ODEProblem(rhs, u0, tspan)
    sol = solve(ode, Tsit5())
    plot(sol)
end



# # Notes for debugging
# function rhs!(du, u, p, t)
#     du = -u
# end

# global tmp[] = (; ode, sol)

# tmp = Ref{Any}()



# ## Digging deeper
#
# If you want to investigate possible bugs deep in a call graph, you
# need to figure out which methods are called under the hood. Julia
# already provides some convenient tools for this in the standard library
# InteractiveUtils (that is loaded by default in a REPL session).

# @which sum([1, 2, 3])
# @edit sum([1, 2, 3])

# You can also use Cthulhu.jl for this purpose:

# using Cthulhu
# @descend sum([1, 2, 3])



# ## Further resources
#
# - [Julia website](https://julialang.org)
# - [Julia documentation](https://docs.julialang.org/en/v1/)
# - [Modern Julia Workflows](https://modernjuliaworkflows.github.io) (in progress)
# - [JuliaNotes.jl](https://m3g.github.io/JuliaNotes.jl/stable/)
# - [Julia Discourse forum](https://discourse.julialang.org)
