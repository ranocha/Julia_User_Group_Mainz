using ForwardDiff, LinearAlgebra
using SparseArrays, SparseDiffTools
using CairoMakie
using Symbolics
include("../../src/helper_functions.jl")
include("perform_timestepping.jl")
include("../../src/adjoint_functions.jl")
include("../../src/addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

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

"""
    Res_adjoint!(F::AbstractArray, u::AbstractArray, Δ, N, BC)
"""
@views function Res_adjoint!(
    R_adjoint::AbstractVector{H},
    Params_array::Vector{<:AbstractArray{H}},
    U_Tuple::NamedTuple,
    Δ::NTuple,
    N::NTuple,
    N_params::Tuple,
    BC::NamedTuple,
    Params::NamedTuple,
    Δt::Float64 
) where {H<:Number}

    dx = Δ      # grid spacing
    nx = N[1]   # grid size

    
    T = U_Tuple.U    # extract 2D fields

    D = Params_array[1]         # extract parameter fields

    # compute heat flux
    qDx = zeros(eltype(D), nx-1)
    qDx = - D .* diff(T) ./ dx
    
    # compute residual
    F_T = zeros(eltype(D), nx)
    F_T[2:end-1] .= (T[2:end-1] .- Params.T_old[2:end-1]) .* inv(Δt) .+ diff(qDx) ./ dx
    
    # boundary conditions
    F_T[1] = T[1] - BC.BC_left
    F_T[end] = T[end] - BC.BC_right

    if eltype(D) == Float64
        # store stuff for visualisation
        Params.qDx .= qDx
    end

    R_adjoint .= [F_T[:];]
    return nothing
end

function Res_adjoint!(
    R_adjoint::AbstractVector{H}, Params_up::AbstractVector{H}, U_Tuple, Δ, N, N_params, BC, Params, Δt::Float64 
) where {H<:Number}
    return Res_adjoint!(R_adjoint, vec_2_vecarray(Params_up, N_params), U_Tuple, Δ, N, N_params, BC, Params, Δt::Float64 )
end

function main(max_timesteps=2, adj_IDX = 2; tmax=500000, verbose=false, maxit=1000, Jac = nothing, colors = nothing)
    
    function Res_closed!(F,U)
        Res!(F, U, Δ, N, BC, Params, Δt)
    end
    
    # Domain
    nx = 500
    lx = 1.0
    Δt = 1e-1
    x = range(-lx / 2, lx / 2, nx + 1)
    xc = (x[2:end] .+ x[1:(end - 1)]) / 2
    dx = lx / (nx - 1)
    
    # Numerics
    N  = (nx,)
    N_params::NTuple = ((nx-1),)
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
    )

    # calculate the adjoint sensitivities
    sensitivities, λ, J_adjoint = adjoint_sensitivities(Jac, U, Δ, N, N_params, BC, Params, Δt, IDX=adj_IDX)

    return U, Jac, colors,time_vec, sensitivities, λ, J_adjoint, Δ, N, N_params, BC, Params, Δt
end

U, Jac, colors, time_vec, sensitivities, λ, J_adjoint, Δ, N, N_params, BC, Params, Δt = main(2, 1:500)
fig1 = Figure()
ax1 = Axis(fig1[1,1], title = "Sensitivities at Nodes $(1:10:500)")
lines!(ax1, sensitivities[1,:],  color=:blue, linewidth = 2, label = "All nodes")
Legend(fig1[1,2], ax1)
display(fig1)

# call the main function, if the jacobian matrix is not passed as an argument, it will be calculated
fig1 = Figure()
ax1 = Axis(fig1[1,1], title = "Sensitivities at Nodes $(1:10:500)")
all_sens = []
for i in 1:50:500
    sens, λ, J_adjoint = adjoint_sensitivities(Jac, U, Δ, N, N_params, BC, Params, Δt, IDX=i)
    lines!(ax1, sens[1,:], colormap = Reverse(:lajolla), color=i, colorrange = (1,500), linewidth = 2, label = "$i")
    push!(all_sens, sens)
end
Legend(fig1[1,2], ax1)
display(fig1)
