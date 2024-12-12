using ForwardDiff, LinearAlgebra
using SparseArrays, SparseDiffTools
using CairoMakie
using Symbolics, Statistics
include("../../src/helper_functions.jl")
include("perform_timestepping.jl")
include("../../src/addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields
include("../../src/adjoint_functions.jl")

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
    dx, dz = Δ      # grid spacing
    T      = U[1]   # extract 2D fields
    nx, nz = N[1]   # grid size

    # compute heat flux
    # in x-direction
    qDx = zeros(eltype(T), (nx-1, nz))
    qDx .= - lin_int(Params.D, 1) .* diff(T, dims=1) ./ dx

    # in z-direction
    qDz = zeros(eltype(T), (nx, nz-1))
    qDz .= - lin_int(Params.D, 2) .* diff(T, dims=2) ./ dz
    
    # compute residual
    F_T = zeros(eltype(T), (nx, nz))
    F_T[2:end-1, 2:end-1] .= (T[2:end-1,2:end-1] .- Params.T_old[2:end-1,2:end-1]) .+ Δt .* (diff(qDx[:,2:end-1], dims=1) ./ dx .+ diff(qDz[2:end-1,:], dims=2) ./ dz .- Params.S[2:end-1,2:end-1])
    
    # boundary conditions
    F_T[end,:] .= T[end,:]  .- BC.BC_right
    F_T[:,end] .= T[:,end]  .- BC.BC_top
    F_T[:,  1] .= T[:,  1]  .- BC.BC_bottom
    F_T[1,  :] .= T[1,  :]  .- BC.BC_left

    if eltype(T) == Float64
        # store stuff for visualisation
        Params.qDx .= qDx
        Params.qDz .= qDz
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

    dx, dz = Δ      # grid spacing
    nx, nz = N[1]      # grid size

    T = U_Tuple.U    # extract 2D fields
    
    D = Params_array[1]         # extract parameter fields

    # compute heat flux
    # in x-direction
    qDx = zeros(eltype(D), (nx-1, nz))
    qDx .= - lin_int(D, 1) .* diff(T, dims=1) ./ dx

    # in z-direction
    qDz = zeros(eltype(D), (nx, nz-1))
    qDz .= - lin_int(D, 2) .* diff(T, dims=2) ./ dz
    
    # compute residual
    F_T = zeros(eltype(D), (nx, nz))
    F_T[2:end-1, 2:end-1] .= (T[2:end-1,2:end-1] .- Params.T_old[2:end-1,2:end-1]) .+ Δt .* (diff(qDx[:,2:end-1], dims=1) ./ dx .+ diff(qDz[2:end-1,:], dims=2) ./ dz .- Params.S[2:end-1,2:end-1])
    
    # boundary conditions
    F_T[end,:] .= T[end,:]  .- BC.BC_right
    F_T[:,end] .= T[:,end]  .- BC.BC_top
    F_T[:,  1] .= T[:,  1]  .- BC.BC_bottom
    F_T[1,  :] .= T[1,  :]  .- BC.BC_left

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


function main(max_timesteps=2, adj_IDX=2; Δt = nothing, tmax=500, verbose=true, maxit=1000, Jac = nothing, colors = nothing)
    
    function Res_closed!(F,U)
        Res!(F, U, Δ, N, BC, Params, Δt)
    end
    
    # Domain
    nx = 100
    nz = 100
    lx = 1.0
    lz = 1.0
    x = range(-lx / 2, lx / 2, nx + 1)
    xc = (x[2:end] .+ x[1:(end - 1)]) / 2
    z = range(-lz / 2, lz / 2, nz + 1)
    zc = (z[2:end] .+ z[1:(end - 1)]) / 2
    dx = lx / (nx - 1)
    dz = lz / (nz - 1)
    
    if isnothing(Δt)
        Δt = 1e-2
    end

    # Numerics
    N  = ((nx,nz),) 
    N_params::NTuple = ((nx, nz),)
    Δ  = (dx,dz)
    
    # initialize solution fields with gaussian
    T    = ones(nx, nz)

    # initialize source term
    Radc = xc.^2 .+ zc'.^2
    S    = ones(nx, nz)
    S[Radc .< 0.005] .= 1.0

    # initialize helper fields
    T_old = copy(T)
    qDx = zeros(nx-1, nz)
    qDz = zeros(nx, nz-1)
    
    # Physical parameters
    D = ones(nx, nz) .* 1e0
    
    # initialize the Arrays and Tuples that will be passed to the residual function
    Params = (; T, T_old, S, qDx, qDz, xc, x, zc, z, D, Radc)
    BC = (; BC_left = 1.0, BC_right = 1.0, BC_top = 1.0, BC_bottom = 1.0)
    U = [T]
    F = vecarray_2_vec(U)

    # initialize residual
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

    return U, Jac, colors, time_vec, sensitivities, λ, J_adjoint, Δ, N, N_params, BC, Params, Δt
end

U, Jac, colors, time_vec, sensitivities, λ, J_adjoint, Δ, N, N_params, BC, Params, Δt = main(2, 1:10000)

fig1 = Figure()
ax1 = Axis(fig1[1,1], title = "Sensitivities at Nodes $(1:10000)")
hm = heatmap!(ax1, reshape(sensitivities, (100,100)), colormap = :berlin, colorrange = (-maximum(abs.(sensitivities)), maximum(abs.(sensitivities))))
Colorbar(fig1[1,2], hm)
display(fig1)



total_sens = zeros(100,100)

for i in 140:171:10000
    fig1 = Figure()
    ax1 = Axis(fig1[1,1], title = "Sensitivities at Node $i")
    sens, λ, J_adjoint = adjoint_sensitivities(Jac, U, Δ, N, N_params, BC, Params, Δt, IDX=i)
    hm = heatmap!(ax1, reshape(sens, (100,100)), colormap = :berlin, colorrange = (-maximum(abs.(sens)), maximum(abs.(sens))))
    Colorbar(fig1[1,2], hm)
    display(fig1)

    total_sens .+= reshape(sens, (100,100))
end

fig1 = Figure()
ax1 = Axis(fig1[1,1], title = "Total Sensitivities")
hm = heatmap!(ax1, total_sens, colormap = :berlin, colorrange = (-maximum(abs.(total_sens)), maximum(abs.(total_sens))))
Colorbar(fig1[1,2], hm)
display(fig1)