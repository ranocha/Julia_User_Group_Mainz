# # Optimizing Julia code

# This session provides an introduction to optimizing Julia code.
# The examples are developed with Julia v1.10.2.

# First, we install all required packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


# ## Basics
#
# The Julia manual contains basic information about performance.
# In particular, you should be familiar with
#
# - https://docs.julialang.org/en/v1/manual/performance-tips/
# - https://docs.julialang.org/en/v1/manual/style-guide/
#
# if you care about performance.
#
# Some additional ressources are
# - https://github.com/albert-de-montserrat/JuliaMainz23
# - https://m3g.github.io/JuliaNotes.jl/stable/
#


# ## Example: A linear advection simulation
#
# This is a classical partial differential equation (PDE) simulation setup.
# If you are not familiar with it, just ignore what's going on - but we need
# an example. This one is similar to several problematic cases I have seen
# in practice.

# First, we generate initial data and store it in a file.
using HDF5
x = range(-1.0, 1.0, length = 1000)
dx = step(x)
h5open("initial_data.h5", "w") do io
    u0 = sinpi.(x)
    write_dataset(io, "x", collect(x))
    write_dataset(io, "u0", u0)
end

# Next, we write our simulation code as a function - as we have learned from
# the performance tips!
function simulate()
    x, u0 = h5open("initial_data.h5", "r") do io
        read_dataset(io, "x"), read_dataset(io, "u0")
    end

    t = 0
    u = u0
    while t < t_end
        u = u + dt * rhs(u)
        t = t + dt
    end

    return t, x, u, u0
end

function rhs(u)
    du = similar(u)
    for i in eachindex(u)
        im1 = i == firstindex(u) ? lastindex(u) : (i - 1)
        du[i] = -(u[i] - u[im1]) / dx
    end
    return du
end

# Now, we can define our parameters, run a simulation,
# and plot the results.
using Plots, LaTeXStrings

t_end = 2.5
dt = 0.9 * dx
t, x, u, u0 = simulate()
plot(x, u0; label = L"u_0", xguide = L"x")
plot!(x, u; label = L"u")

# Maybe you can already spot problems in the code above.
# Anyway, before you begin optimizing some code, you should
# test it and make sure it's correct. Let's assume this is the
# case here. You should write tests making sure that the code
# continues to work as expected while we optimize it. We will
# not do this here right now.


# ## Profiling
#
# First, we need to measure the performance of our code to
# see whether it's already reasonable. For this, you can use
# `@time` - or better BenchmarkTools.jl.

using BenchmarkTools
@benchmark simulate()

# This shows quite a lot of allocations - typically a bad sign
# if you are working with numerical simulations.

# There are also profiling tools available in Julia.
# Some of them are already loaded in the Visual Studio Code
# extension for Julia. If you prefer working in the REPL,
# you can install ProfileView.jl (`@profview simulate()`) and
# PProf.jl (`pprof()` after creating a profile via `@profview`).

@profview simulate()

# Task: Optimize the code!





# ## Introduction to generic Julia code and AD
#
# One of the strengths of Julia is that you can use quite a few
# modern tools like AD. However, you need to learn Julia a bit
# to use everything efficiently.
# Here, we concentrate on ForwardDiff.jl.

using ForwardDiff
using Statistics

function foo0(x)
    vector = zeros(typeof(x), 100)
    foo0!(vector, x)
end

function foo0!(vector, scalar)
    for i in eachindex(vector)
        vector[i] = atan(i * scalar) / π
    end

    for _ = 1:1000
        for i in eachindex(vector)
            vector[i] = cos(vector[i])
        end
    end

    return mean(vector)
end

let x = 2 * π
    @show foo0(x)
    @show ForwardDiff.derivative(foo0, x)
    @benchmark foo0($x)
end


# The code above is basically a fixed point iteration.
# By doing some analysis, one could figure out that it
# should converge to the solution `x` of `x == cos(x)`,
# the "Dottie number". See https://en.wikipedia.org/wiki/Dottie_number

using SimpleNonlinearSolve

function dottie_number()
    f(u, parameters) = cos(u) - u
    bounds = (0.0, 2.0)
    prob = IntervalNonlinearProblem(f, bounds)
    sol = solve(prob, ITP())
    return sol.u
end

dottie_number()


# Next, we introduce a custom struct. Can you spot the problems?

struct MyData1
    scalar
    vector
end

function foo1(x)
    data = MyData1(x, zeros(typeof(x), 100))
    foo1!(data)
end

function foo1!(data)
    (; scalar, vector) = data

    for i in eachindex(vector)
        vector[i] = atan(i * scalar) / π
    end

    for _ = 1:1000
        for i in eachindex(vector)
            vector[i] = cos(vector[i])
        end
    end

    return mean(vector)
end

let x = 2 * π
    @show foo1(x)
    @show ForwardDiff.derivative(foo1, x)
    @benchmark foo1($x)
end


# We can fix type instabilities by introducing types explicitly.
# But we lose the possibility to use AD this way!

struct MyData2
    scalar::Float64
    vector::Vector{Float64}
end

function foo2(x)
    data = MyData2(x, zeros(typeof(x), 100))
    foo1!(data)
end

let x = 2 * π
    @show foo2(x)
    @benchmark foo2($x)
end

let x = 2 * π
    ForwardDiff.derivative(foo2, x)
end

# Thus, we need parametric types!

struct MyData3{T}
    scalar::T
    vector::Vector{T}
end

function foo3(x)
    data = MyData3(x, zeros(typeof(x), 100))
    foo1!(data)
end

let x = 2 * π
    @show foo3(x)
    @show ForwardDiff.derivative(foo3, x)
    @benchmark foo3($x)
end


# We could also use something like

struct MyData4{T1 <: Real, T2 <: AbstractVector{T1}}
    scalar::T1
    vector::T2
end


function foo(x, y)
    # ....
    return x + sum(y)
end

function bar(x::Real, y::AbstractVector{<:Real})
    # ....
    return x + sum(y)
end
