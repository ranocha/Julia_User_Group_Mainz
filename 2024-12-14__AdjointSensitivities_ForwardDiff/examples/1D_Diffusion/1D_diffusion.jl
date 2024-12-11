using ForwardDiff, LinearAlgebra
using SparseArrays, SparseDiffTools
using CairoMakie
using Symbolics
include("/Users/jacob/Documents/GitHub/AD_FD_examples/src/helper_functions.jl")
include("/Users/jacob/Documents/GitHub/AD_FD_examples/examples/1D_Diffusion/perform_timestepping.jl")
include("/Users/jacob/Documents/GitHub/AD_FD_examples/src/addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

"""
    Res!(F::AbstractArray, u::AbstractArray, Δ, N, BC)
"""
@views function Res!(
    Fvec::AbstractVector{H},
    U::Vector{<:AbstractArray{H}},
    Δ::NTuple,
    N::NTuple,
    BC::NamedTuple,
    Params::NamedTuple,
    Δt,
) where {H<:Number}
    dx = Δ      # grid spacing
    T  = U[1]   # extract 2D fields
    nx = N[1]   # grid size

    # compute heat flux
    qDx = zeros(eltype(T), nx-1)
    qDx = - Params.D .* diff(T) ./ dx
    
    # compute residual
    F_T = zeros(eltype(T), nx)
    F_T[2:end-1] .= (T[2:end-1] .- Params.T_old[2:end-1]) .* inv(Δt) .+ diff(qDx) ./ dx
    
    # boundary conditions
    F_T[1] = T[1] - BC.BC_left
    F_T[end] = T[end] - BC.BC_right

    if eltype(T) == Float64
        # store stuff for visualisation
        Params.qDx .= qDx
    end

    Fvec .= [F_T[:];]
    return nothing
end

function Res!(
    Fup::AbstractVector{H}, up::AbstractVector{H}, Δ, N, BC, Params, Δt
) where {H<:Number}
    return Res!(Fup, vec_2_vecarray(up, N), Δ, N, BC, Params, Δt)
end

function main(max_timesteps=2; tmax=500000, verbose=false, maxit=1000, Jac = nothing, colors = nothing)
    
    function Res_closed!(F,U)
        Res!(F, U, Δ, N, BC, Params, Δt)
    end
    
    # Domain
    nx = 500
    lx = 1.0
    Δt = 1e-2
    x = range(-lx / 2, lx / 2, nx + 1)
    xc = (x[2:end] .+ x[1:(end - 1)]) / 2
    dx = lx / (nx - 1)
    
    # Numerics
    N = (nx,) 
    Δ  = (dx,)
    
    # initialize solution fields
    T   = ones(nx)
    T[1] = 5.0
    
    # initialize helper fields
    T_old = copy(T)
    qDx = zeros(nx-1)
    
    # Physical parameters
    k   = ones(nx)   .* 3.0
    ρ   = ones(nx)   .* 3000
    cp  = ones(nx)   .* 1000
    D    = ones(nx-1) .* 1e0
    
    # initialize the Arrays and Tuples that will be passed to the residual function
    Params = (;k, ρ, cp, T_old, qDx, xc, x, D)
    BC = (; BC_left = 5.0, BC_right = 1.0)
    U = [T]
    F = vecarray_2_vec(U)

    Res_closed!(F, U)

    if isnothing(Jac)
        println("--------------------------------\nCalculating the Jacobian matrix...\n")

        # Estimate the sparsity pattern using a brute force approach
        Uvec        =   vecarray_2_vec(U)
        J2          =   ForwardDiff.jacobian(Res_closed!, F, Uvec); # slow and non-sparse jacobian; only needs to be done once
        Jac         =   sparse(Float64.(abs.(J2).>0))
        colors      =   matrix_colors(Jac) 

        println("--------------------------------\nJacobian matrix calculated.\n")
    end

    time_vec, U = perform_timestepping(
        Params, 
        U, 
        F, 
        Jac, 
        colors, 
        Δt, 
        Res_closed!, 
        N;
        tmax=tmax,
        max_timesteps=max_timesteps,
        verbose=verbose,
        maxit=maxit,
        save_fig=true,
        fig_name="1D_diffusion"
    )

    return U, Jac, colors
end


# call the main function, if the jacobian matrix is not passed as an argument, it will be calculated

U, Jac, colors = main()
