# example_stokes.jl
#
# We solve the 2D stokes equations using a staggered grid formulation with Vx,Vz,P as unknowns
# The equations are:
#        ∂Vx/∂x + ∂Vz/∂z = -P/K
#   -∂P/∂x + ∂τxx/∂x + ∂τxz/∂x = 0
#   -∂P/∂z + ∂τxz/∂x + ∂τzz/∂x = ρg
#
#   τxx = 2 η εxx 
#   τzz = 2 η εzz 
#   τxz = 2 η εxz 
#
#   εxx = ∂Vx/∂x
#   εzz = ∂Vz/∂z
#   εxz = (∂Vx/∂z + ∂Vz/∂x)/2

using SparseDiffTools, LinearAlgebra, GeoParams, Statistics
using Symbolics
using GLMakie, JLD2
using GeophysicalModelGenerator
Makie.inline!(false)

include("./src/addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

update_density!(ρ, phase, rheology) = update_field!(ρ, phase, rheology, compute_density)
update_viscosity!(η, phase, rheology) = update_field!(η, phase, rheology, compute_viscosity)

function update_field!(A, phase, rheology, f::F) where F
    args = (;)
    for i in eachindex(A)
        A[i] = f(rheology, phase[i], args)
    end
    return nothing
end

average_2D(x) = (x[2:end,2:end] + x[1:end-1,2:end] + x[2:end,1:end-1] + x[1:end-1,1:end-1])/4

function lin_int(A, dim)
    if dim == 1
        A = (A[1:(end - 1), :] .+ A[2:end, :]) ./ 2.0
    elseif dim == 2
        A = (A[:, 1:(end - 1)] .+ A[:, 2:end]) ./ 2.0
    end
    return A
end

function harmonic_avg(η1, η2)
    return 2.0 * (η1 * η2) / (η1 + η2)
end

function smooth_viscosity_harmonic(η, Nx, Nz)
    for i in 2:Nx-1
        for j in 2:Nz-1
            η[i, j] = harmonic_avg(
                harmonic_avg(η[i, j], η[i+1, j]),
                harmonic_avg(η[i, j], η[i, j+1])
            )
        end
    end
    return η
end

function smooth_viscosity(η, Nx, Nz)
    η_smooth = copy(η)
    for i in 2:Nx-1
        for j in 2:Nz-1
            η_smooth[i, j] = 0.2 * (η[i, j] + η[i+1, j] + η[i-1, j] + η[i, j+1] + η[i, j-1])
        end
    end
    η_smooth[:,1] .= η_smooth[:,2]
    η_smooth[:,end] .= η_smooth[:,end-1]
    η_smooth[1,:] .= η_smooth[2,:]
    η_smooth[end,:] .= η_smooth[end-1,:]
    return η_smooth
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function update_stresses!(τxx, τzz, τxz, εxx,εzz,εxz, P, Params, Δt)
    I1   = CartesianIndex(1,0)
    I2   = CartesianIndex(1,1)
    I3   = CartesianIndex(0,1)
    τxz .= 0.0
    for I in CartesianIndices(εxx)
        # old stresses
        τ_o = (
            Params.τxx_old[I],
            Params.τzz_old[I],
            (Params.τxz_old[I], Params.τxz_old[I+I1], Params.τxz_old[I+I2], Params.τxz_old[I+I3]),
        )
        ε_phase = (εxx[I], εzz[I], (εxz[I], εxz[I+I1], εxz[I+I2], εxz[I+I3]))
        phase_i = (Params.phasec[I], Params.phasec[I], (Params.phasev[I], Params.phasev[I+I1], Params.phasev[I+I2], Params.phasev[I+I3]))
        
        args = (dt=Δt, P = P[I])
        τij, τII = compute_τij(rheology, ε_phase, args, τ_o, phase_i)
        
        τxx[I   ]  = τij[1]
        τzz[I   ]  = τij[2]
        # we average shear stresses @ vertices
        τxz[I   ] += τij[3]*0.25;
        τxz[I+I1] += τij[3]*0.25;
        τxz[I+I2] += τij[3]*0.25;
        τxz[I+I3] += τij[3]*0.25;
    end

    return nothing
end

function update_stresses_visc!(
    τxx, τzz, τxz, 
    εxx, εzz, εxz, 
    ηv,
    P
)
    Nx, Nz = size(εxx)

    # Reset shear stresses at the vertices before updating
    τxz .= 0.0  

    ηc = average_2D(ηv)
    # Update normal stresses τxx and τzz
    for i in 1:Nx
        for j in 1:Nz
            local_viscosity = ηc[i,j]
            τxx[i,j] = 2.0 * local_viscosity * εxx[i,j] - (1.0 / 3.0) * P[i,j]
            τzz[i,j] = 2.0 * local_viscosity * εzz[i,j] - (1.0 / 3.0) * P[i,j]
        end
    end

    # Update shear stress τxz at the vertices with averaged viscosity
    for i in 1:Nx
        for j in 1:Nz
            local_viscosity_avg = 0.25 * (ηv[i, j] + ηv[i+1, j] + ηv[i, j+1] + ηv[i+1, j+1])
            τxz[i, j]       += 2.0 * local_viscosity_avg * εxz[i,j] * 0.25
            τxz[i+1, j]     += 2.0 * local_viscosity_avg * εxz[i,j] * 0.25
            τxz[i, j+1]     += 2.0 * local_viscosity_avg * εxz[i,j] * 0.25
            τxz[i+1, j+1]   += 2.0 * local_viscosity_avg * εxz[i,j] * 0.25
        end
    end

    # Handle the boundary conditions explicitly for shear stress (free-slip BC)
    τxz[1, :] .= 0.0
    τxz[:, 1] .= 0.0
    τxz[end, :] .= 0.0
    τxz[:, end] .= 0.0

    return nothing
end



# START Residual function for the forward problem -----------------------------------

"""
    Res!(F::AbstractArray, u::AbstractArray, Δ, N, BC)
"""
@views function Res!(Fvec::AbstractVector{T}, U::Vector{<:AbstractArray{T}}, Δ::NTuple, N::NTuple, BC::NamedTuple, Params::NamedTuple, rheology, Δt) where T<:Number

    dx, dz      = Δ    # grid spacing
    Vx, Vz, P   = U    # extract 2D fields
    Nx, Nz      = N[3] # grid size
    
    Plocal = copy(P) .- Params.P_shift[1] # set average P to zero @ top

    # Strainrates
    εvol = diff(Vx,dims=1)./dx + diff(Vz,dims=2)./dz
    εxx  = diff(Vx,dims=1)./dx - 1/3*εvol
    εzz  = diff(Vz,dims=2)./dz - 1/3*εvol
    εxz  = zeros(eltype(Vx), Nx+1, Nz+1)
    εxz[2:end-1,2:end-1] = (diff(Vx[2:end-1,:],dims=2)./dz + diff(Vz[:,2:end-1],dims=1)./dx)
    
    # Update deviatoric stresses (using GeoParams)
    τxx = zeros(eltype(Vx), Nx, Nz)
    τzz = zeros(eltype(Vx), Nx, Nz)
    τxz = zeros(eltype(Vx), Nx+1, Nz+1)

    # update stresses usinf GeoParams
    # update_stresses!(τxx, τzz, τxz, εxx,εzz,εxz, P, Params, Δt)

    # by hand update stresses
    update_stresses_visc!(τxx, τzz, τxz, εxx, εzz, εxz, Params.ηv, Plocal)

    # Free slip BC's
    τxz[1,  :] .= 0
    τxz[:,end] .= 0
    τxz[:,  1] .= 0
    τxz[end,:] .= 0

    update_density!(Params.ρ,    Params.phasec, rheology)
  
    # Mass balance
    Fmass = zero(P)
    Fmass = εxx +  εzz #+ Params.K.*P 

    # x-momentum balance
    Fforce1             = zero(Vx)
    Fforce1[2:end-1,:]  = diff(-Plocal + τxx, dims=1)./dx + diff(τxz[2:end-1,:], dims=2)./dz
    Fforce1[1,:]        = Vx[1,:]       #.- BC.εxx_bg*Params.x[1];
    Fforce1[end,:]      = Vx[end,:]     #.- BC.εxx_bg*Params.x[end];

    # z-momentum balance
    Fforce2             =   zero(Vz)
    g                   =   rheology[1].Gravity[1].g.val;
    Fforce2[:,2:end-1]  =   diff(τxz[:,2:end-1], dims=1)./dx + diff(-Plocal + τzz, dims=2)./dz - g*lin_int(ρ, 2) 
    Fforce2[:,      1]  =   Vz[:,1]     #.- BC.εzz_bg*Params.z[1];
    Fforce2[:,    end]  =   Vz[:,end]   #.- BC.εzz_bg*Params.z[end];
    
    if eltype(Vx) == Float64
        τxz2_c = average_2D(τxz.^2)
        εxz2_c = average_2D(εxz.^2)
        
        Params.τxx .= τxx
        Params.τzz .= τzz
        Params.τxz .= τxz
        Params.εII .= sqrt.(0.5*εxx.^2 + 0.5.*εzz.^2 + εxz2_c)
        Params.τII .= sqrt.(0.5*τxx.^2 + 0.5.*τzz.^2 + τxz2_c)
    end

    Fvec .= [Fforce1[:]; Fforce2[:]; Fmass[:]]
    return nothing
end

function Res!(Fup::AbstractVector{T}, up::AbstractVector{T}, Δ, N, BC, Params, rheology, Δt) where {T<:Number}
    return Res!(Fup, vec_2_vecarray(up, N), Δ, N, BC, Params, rheology, Δt)
end

Res_closed! = (F,U) -> Res!(F, U, Δ, N, BC, Params, rheology, Δt)        # create a function with only 1 input parameter

# END Residual function for the forward problem -------------------------------------
# START Residual function for the adjoint problem -----------------------------------

@views function Res_adjoint!(
    R_adjoint::AbstractVector{T},
    Params_array::Vector{<:AbstractArray{T}},
    U_Tuple::NamedTuple,
    Δ::NTuple,
    N::NTuple,
    BC::NamedTuple,
    Params::NamedTuple                                                                 
) where {T<:Number}

   
    dx, dz      = Δ    # grid spacing
    Nx, Nz      = N[3] # grid size
    Vx, Vz, P   = U_Tuple.U    # extract 2D fields
    
    # extract parameters we would like to differentiate for
    ρ, η = Params_array
    
    Plocal = copy(P) .- Params.P_shift[1] # set average P to zero @ top

    # Strainrates
    εvol = diff(Vx,dims=1)./dx + diff(Vz,dims=2)./dz
    εxx  = diff(Vx,dims=1)./dx - 1/3*εvol
    εzz  = diff(Vz,dims=2)./dz - 1/3*εvol
    εxz  = zeros(eltype(Vx), Nx+1, Nz+1)
    εxz[2:end-1,2:end-1] = (diff(Vx[2:end-1,:],dims=2)./dz + diff(Vz[:,2:end-1],dims=1)./dx)

    # Update deviatoric stresses (using GeoParams)
    τxx = zeros(eltype(ρ), Nx, Nz)
    τzz = zeros(eltype(ρ), Nx, Nz)
    τxz = zeros(eltype(ρ), Nx+1, Nz+1)

    # update stresses usinf GeoParams
    # update_stresses!(τxx, τzz, τxz, εxx,εzz,εxz, P, Params, Δt)

    # by hand update stresses
    update_stresses_visc!(τxx, τzz, τxz, εxx, εzz, εxz, η, Plocal)

    # Free slip BC's
    τxz[1,  :] .= 0
    τxz[:,end] .= 0
    τxz[:,  1] .= 0
    τxz[end,:] .= 0

    # update_density!(ρ,    Params.phasec, rheology)

    # Mass balance
    Fmass = zeros(eltype(ρ), size(P))
    Fmass = εxx +  εzz #+ Params.K.*P 

    # x-momentum balance
    Fforce1             = zeros(eltype(ρ), size(Vx))
    Fforce1[2:end-1,:]  = diff(-Plocal + τxx, dims=1)./dx + diff(τxz[2:end-1,:], dims=2)./dz
    Fforce1[1,:]        = Vx[1,:]       #.- BC.εxx_bg*Params.x[1];
    Fforce1[end,:]      = Vx[end,:]     #.- BC.εxx_bg*Params.x[end];

    # z-momentum balance
    Fforce2             =   zeros(eltype(ρ), size(Vz))
    g                   =   rheology[1].Gravity[1].g.val;
    Fforce2[:,2:end-1]  =   diff(τxz[:,2:end-1], dims=1)./dx + diff(-Plocal + τzz, dims=2)./dz - g*lin_int(ρ, 2)  
    Fforce2[:,      1]  =   Vz[:,1]     #.- BC.εzz_bg*Params.z[1];
    Fforce2[:,    end]  =   Vz[:,end]   #.- BC.εzz_bg*Params.z[end];

    if eltype(ρ) == Float64
        τxz2_c = average_2D(τxz.^2)
        εxz2_c = average_2D(εxz.^2)
        
        Params.τxx .= τxx
        Params.τzz .= τzz
        Params.τxz .= τxz
        Params.εII .= sqrt.(0.5*εxx.^2 + 0.5.*εzz.^2 + εxz2_c)
        Params.τII .= sqrt.(0.5*τxx.^2 + 0.5.*τzz.^2 + τxz2_c)
    end

    R_adjoint .= [Fforce1[:]; Fforce2[:]; Fmass[:]]
    return nothing
end

function Res_adjoint!(
    R_adjoint::AbstractVector{T}, Params_up::AbstractVector{T}, U_Tuple, Δ, N, BC, Params
) where {T<:Number}
    return Res_adjoint!(R_adjoint, vec_2_vecarray(Params_up, N_params), U_Tuple, Δ, N, BC, Params)
end

Res_adjoint_closed! = (R_adjoint, Params_array) -> Res_adjoint!(R_adjoint, Params_array, U_Tuple, Δ, N, BC, Params)

function compute_λ(J_forward, ∇_F)
    return - J_forward' \ ∇_F
end

function dR_dp(Params_array, U)
    # Ensure R_adjoint is correctly sized based on the forward problem residual
    R_adjoint = similar(vecarray_2_vec(U))
    
    # Calculate the adjoint residual for the given Params_array
    Res_adjoint_closed!(R_adjoint, Params_array)
    
    # Compute the Jacobian of the adjoint residual with respect to the parameters
    J_adjoint = ForwardDiff.jacobian(Res_adjoint_closed!, R_adjoint, vecarray_2_vec(Params_array))
    
    return J_adjoint
end

function adjoint_sensitivities(J_forward, U; IDX=100, IDX1=nothing)

    println("Computing sensitivities for node(s) $IDX")

    global U_Tuple = (U=U,)
    Params_array = [Params.ρ, Params.ηv]
    R_adjoint = vecarray_2_vec(U)

    # Compute the gradient of the objective function w.r.t. U
    # in the case that we want the sensitivites this is just a mask with 1's and 0's

    ∇_F = zeros(size(R_adjoint))
    if typeof(IDX) == Int
        ∇_F[IDX] = 1.0  # Sensitivity with respect to the node #
    else
        ∇_F[IDX] .= 1.0  # Sensitivity with respect to the node #
    end

    if IDX1 != nothing
        ∇_F[IDX1] .= 1.0  # Sensitivity with respect to the node #
    end

    # Solve the adjoint problem J^T * lambda = -grad_F
    λ = compute_λ(J_forward, ∇_F)

    # println("λ: ", λ)

    # Compute dR/dp (Jacobian of residual w.r.t. parameters)
    J_adjoint = dR_dp(Params_array, U)

    # println("J_adjoint: ", J_adjoint)

    # Compute the sensitivity dF/dp
    sensitivities = -λ' * J_adjoint

    #normnalize sensitivities
    # sensitivities = sensitivities ./ norm(sensitivities)
    # reshape sensitivities to 2D array
    sensitivities_ρ = reshape(sensitivities[1:length(U[3])],     size(Params_array[1]))
    sensitivities_η = reshape(sensitivities[length(U[3])+1:end], size(Params_array[2]))

    sensitivities = (sensitivities_ρ, sensitivities_η)

    return sensitivities, λ, J_adjoint
end

# END Residual function for the adjoint problem -----------------------------------

function perform_timestepping(Params, U, F, Jac, colors, Δt;  tmax=1e3, max_timesteps=1e6, verbose=false, maxit=1000)
    time = 0.0
    it = 0
    global Δt

    τxx_time = Float64[]
    time_vec = Float64[]
    time     = 0;
    Usol     = copy(U)
    U0       = copy(U)

    while time<=tmax && it<max_timesteps
        converged = false
        maxit = 100;
        copyto!(U0, U)
        while !converged
            Usol, converged, its = nonlinear_solution(F, U, Jac, colors, maxit=maxit, verbose=verbose);

            if !converged
                copyto!(U, U0)
                Δt0 = Δt
                Δt  = Δt / 2
                println("  Did NOT converged in $its steps; decreasing Δt from $Δt0 to $(Δt)")
                
            elseif its<10  #&& Δt<max_Δt
                println("  Converged in $its steps; Δt=$Δt")
                Δt *= 1.2
            end
        end
        time += Δt
        it   += 1
        

        # plot surface velocities
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel = "x [km]", ylabel = "z [km]", title = "Surface Vz")
        lines!(ax, Params.xc*CharUnits.length.val, Usol[2][:,end-32]*CharUnits.length.val, color = :blue)
        display(fig)

        # Update old parameters
        copyto!(U, Usol)
        Params.τxx_old   .= Params.τxx
        Params.τzz_old   .= Params.τzz
        Params.τxz_old   .= Params.τxz
        Params.P_shift[1] = mean(U[3][:,end])

        τxx = mean(Params.τxx)*CharUnits.stress.val
        time_dim = time*CharUnits.time.val
        push!(τxx_time, τxx)
        push!(time_vec, time_dim)

        τII_dim = round(mean(Params.τII)*CharUnits.stress.val,digits=2)
        τII_max_dim = round(maximum(Params.τII*CharUnits.stress.val), digits=2)
        
        println("Timestep $it/$max_timesteps, time=$time, Δt=$Δt, τII=$(τII_dim) MPa, max(τII)=$τII_max_dim MPa")
    end
    return time_vec, τxx_time, U
end

# Setup problem
Nx, Nz = 120,60
# parameters in dimensional units
# L_d     =    10km;
# zbot_d  =    -10km
# W_d     =    20km;
# R_d     =    0.25km;
# Δt_d    =    100*yr     # timestep
# tmax_d  =    2Myrs;                         # in Myrs
# K       =    1e-6
# εxx_d   =    1e-100/s                   #1e-15/s;
# εzz_d   =    1e-100/s                   #-1e-15/s;

L_d     =    1.45m
zbot_d  =    -1m
W_d     =    2m
R_d     =    0.25m
Δt_d    =    1*s     # timestep
tmax_d  =    200*s;                         # in Myrs
εxx_d   =    1e-150/s                  #1e-15/s;
εzz_d   =    -1e-150/s                  #-1e-15/s;
K = 1e-6

# CharUnits = GEO_units(length=10km, temperature=1000C, stress=10MPa, viscosity=1e23Pas)
CharUnits = SI_units(length=1m, temperature=1K, stress=1Pa, viscosity=1Pas)
L,W,R,Δt,tmax,εxx_bg,εzz_bg,zbot = nondimensionalize(L_d, W_d, R_d, Δt_d, tmax_d, εxx_d, εzz_d,zbot_d, CharUnits) 

# Numerics
N   =   ((Nx+1,Nz), (Nx,Nz+1), (Nx,Nz))
N_params = ((Nx, Nz), (Nx+1, Nz+1))
Vx  =   rand(N[1]...)
Vz  =   rand(N[2]...)
P   =   rand(N[3]...)
x   =   range(-W/2,W/2,Nx+1)
z   =   range(zbot  ,zbot+L  ,Nz+1)
xc  =   (x[2:end] .+ x[1:end-1])/2
zc  =   (z[2:end] .+ z[1:end-1])/2

# Material fields
ρ   =    ones(Nx,Nz) .* 1.2             # density
η   =    ones(Nx,Nz)             # viscosity @ center
ηv  =    ones(Nx+1,Nz+1)         # viscosity @ vertices
K   =    ones(Nx, Nz) .* 1e-6    # bulk modulus

# Helper fields
εxx     =  zeros(Nx,Nz)
εzz     =  zeros(Nx,Nz)
εII     =  zeros(Nx,Nz)
εxz     =  zeros(Nx+1,Nz+1)
τxx     =  zeros(Nx,Nz)
τzz     =  zeros(Nx,Nz)
τII     =  zeros(Nx,Nz)
τxz     =  zeros(Nx+1,Nz+1)
τxx_old =  zeros(Nx,Nz)
τzz_old =  zeros(Nx,Nz)
τxz_old =  zeros(Nx+1,Nz+1)
phasec  =  ones(Int64, Nx,Nz)      # material phase tag @ center
phasev  =  ones(Int64, Nx+1,Nz+1)  # material phase tag @ vertices

dx  = x.step.hi
dz  = z.step.hi
Δ   = (dx,dz)

# Set material phases
# Define the dimensions of the box in phase 2

XC, ZC = meshgrid(xc, zc)
XV, ZV = meshgrid(x, z)

# defien the subduction zone
# Set material phases
# Define the dimensions of the box in phase 2
box_width = 0.1  # Half-width of the box (in non-dimensional units)
box_height = 0.1 # Half-height of the box (in non-dimensional units)

# Iterate over the center points
for (i, xv) in enumerate(xc), (j, zv) in enumerate(zc)
    # Define the box in the middle of the domain with phase 2
    if abs(xv) < box_width && abs(zv - (zbot + (L-0.5)/2)) < box_height
        phasec[i, j] = 2
        ρ[i, j] = 1.5
    else
        phasec[i, j] = 1
        ρ[i, j] = 1
    end
end

# Iterate over the vertex points
for (i, xv) in enumerate(x), (j, zv) in enumerate(z)
    # Define the box in the middle of the domain with phase 2
    if abs(xv) < box_width && abs(zv - (zbot + (L-0.5)/2)) < box_height
        phasev[i, j] = 2
        ηv[i, j] = 1e3
    else
        phasev[i, j] = 1
        ηv[i, j] = 1
    end
end

ηv[:,end-31] .= 0.9
ηv[:,end-30] .= 0.8
ηv[:,end-29] .= 0.7
ηv[:,end-28] .= 0.6
ηv[:,end-27] .= 0.5
ηv[:,end-26] .= 0.4
ηv[:,end-25] .= 0.3
ηv[:,end-24] .= 0.2
ηv[:,end-23:end] .= 0.1

# smooth the viscosity
# η = smooth_viscosity_harmonic(η, Nx, Nz)
# ηv = smooth_viscosity_harmonic(ηv, Nx+1, Nz+1)
ηv = smooth_viscosity(ηv, Nx+1, Nz+1)
ηv = smooth_viscosity(ηv, Nx+1, Nz+1)
K  = smooth_viscosity(K, Nx, Nz)
K  = smooth_viscosity(K, Nx, Nz)
# GeoParams rheology
rheology = (
        SetMaterialParams(;
            Name              = "Matrix",
            Gravity             = ConstantGravity(g=1.0m/s^2),
            Phase             = 1,
            Density           = ConstantDensity(ρ=1.0kg/m^3),
            CompositeRheology = CompositeRheology((LinearViscous(η=1e0Pas),
                                                   #ConstantElasticity(),
                                                   #DruckerPrager_regularised(ϕ=30, C=80MPa, η_vp=1e19Pas)
                                                   )),
            CharDim             = CharUnits    
        ),
        SetMaterialParams(;
            Name              = "Inclusion",
            Gravity             = ConstantGravity(g=1.0m/s^2),
            Phase             = 2,
            Density           = ConstantDensity(ρ=1.5kg/m^3),
            CompositeRheology = CompositeRheology((LinearViscous(η=1e3Pas),
                                                   #ConstantElasticity(),
                                                   #DruckerPrager_regularised(ϕ=30, C=80MPa, η_vp=1e19Pas)
                                                   )),
            CharDim             = CharUnits    
        ),
        SetMaterialParams(;
            Name              = "StickyAir",
            Gravity             = ConstantGravity(g=1.0m/s^2),
            Phase             = 3,
            Density           = ConstantDensity(ρ=1000kg/m^3),
            
            CompositeRheology = CompositeRheology((LinearViscous(η=1e18Pas),
                                                    ConstantElasticity(),
                                                    #DruckerPrager_regularised(ϕ=30, C=80MPa, η_vp=1e19Pas)
                                              )),
            CharDim           = CharUnits,    
            ),
)

P_shift = [0.0]
Params  = (; K, εxx,εzz,εxz,τxx,τzz,τxz,τxx_old,τzz_old,τxz_old, ρ,  phasec, phasev, xc,zc,x,z, Δt, εII, τII, P_shift, ηv)                 # additional parameters of the PDE
BC      = (; εxx_bg,εzz_bg)
U       = [Vx,Vz,P]
F       = vecarray_2_vec(U)


Res_closed!(F,U)

# Sparsity pattern of jacobian
#@time sparsity    =   jacobian_sparsity(Res_closed!,F,U)
#J           =   Float64.(sparsity)
#colors      =   matrix_colors(J) 

# Estimate the sparsity pattern using a brute force approach
using ForwardDiff, SparseArrays
Uvec        =   vecarray_2_vec(U)
J2          =   ForwardDiff.jacobian(Res_closed!, F, Uvec); # slow and non-sparse jacobian; only needs to be done once
Jac         =   sparse(Float64.(abs.(J2).>0))
colors      =   matrix_colors(Jac) 

# initial guesses for Vx,Vz,P
for (i,xv) in enumerate(x), (j,zv) in enumerate(zc)
    Vx[i,j] = BC.εxx_bg*xv
end
for (i,xv) in enumerate(xc), (j,zv) in enumerate(z)
    Vz[i,j] = BC.εzz_bg*zv
end
for (i,xv) in enumerate(xc), (j,zv) in enumerate(zc)
    P[i,j] = rheology[1].Density[1].ρ.val* rheology[1].Gravity[1].g.val
end

U   = [Vx,Vz,P]

if 1==0
    #jldsave("U.jld2"; U, τxz_old=Params.τxz_old, τxx_old=Params.τxx_old, τzz_old=Params.τzz_old);
    U,τxz_old,τxx_old,τzz_old = load("U.jld2", "U", "τxz_old","τxx_old","τzz_old")
    Params.τxx_old .= τxx_old
    Params.τzz_old .= τzz_old
    Params.τxz_old .= τxz_old
end

time_vec, τxx_time, Usol = perform_timestepping(Params, U, F, Jac, colors, Δt;  
                        max_timesteps=1, tmax=tmax, verbose=true)


# plot
#lines(time_vec, τxx_time, color=:blue, axis=(xlabel="Time [Myrs]", ylabel="τxx [MPa]"))

data = Params.τII*CharUnits.stress.val; title_str="τII [MPa]"

P = (Usol[3] .- Params.P_shift[1])*CharUnits.stress.val
#data = P ; title_str="P [MPa]"

fig, ax, hm = heatmap(xc*CharUnits.length.val,zc*CharUnits.length.val,data)
Colorbar(fig[1,2], limits = extrema(data), label=title_str)
fig


sensitivities, λ, J_adjoint = adjoint_sensitivities(J2, Usol, IDX=Nx*Nz+Nx*(Nz-32)-Int(7*Nx/10):Nx*Nz+Nx*(Nz-32)-Int(3*Nx/10))

IDX_array = [i for i in 1:1:length(Usol[3])]
IDX_array = reshape(IDX_array, size(Usol[3]))
IDX = []
for i in Nx*(Nz-20)-Int(7*Nx/10):Nx*(Nz-20)-Int(3*Nx/10)
    IDX_plot = findall(IDX_array .== i)
    push!(IDX, IDX_plot)
end


#plot same sections as Georg
ind_box = findall(XC .> 0.9 .|| XC .< -0.9 .|| ZC .< -0.8 .|| ZC .> -0.2)
sens_ρ = copy(sensitivities[1])
sens_η = copy(sensitivities[2])
for i in eachindex(ind_box)
    sens_ρ[ind_box[i]] = 0
    sens_η[ind_box[i]] = 0
end


fig = Figure(size=(1000, 800))
#Vz
ax = Axis(fig[1, 1], xlabel = "x [km]", ylabel = "z [km]", title = "Vz")
hm = heatmap!(ax, xc*CharUnits.length.val, zc[1:end-32]*CharUnits.length.val, Usol[2][:,1:end-32]*CharUnits.velocity.val, colormap = :viridis)
# arrows!(ax, xc*CharUnits.length.val, zc*CharUnits.length.val, lin_int(Usol[1],1), lin_int(Usol[2],2), lengthscale = 0.005, arrowsize = 1)
Colorbar(fig[1, 2], hm, label = "Vz [m/s]")
# #Vx
ax2 = Axis(fig[1, 3], xlabel = "x [km]", ylabel = "z [km]", title = "Vx")
hm2 = heatmap!(ax2, xc*CharUnits.length.val, zc[1:end-32]*CharUnits.length.val, Usol[1][:,1:end-32]*CharUnits.velocity.val, colormap = :viridis)
# arrows!(ax2, xc*CharUnits.length.val, zc*CharUnits.length.val, lin_int(Usol[1],1), lin_int(Usol[2],2), lengthscale = 0.005, arrowsize = 1)
Colorbar(fig[1, 4], hm2, label = "Vx [m/s]")
#P
ax3 = Axis(fig[2,1], xlabel = "x [km]", ylabel = "z [km]", title = "Stress")
hm3 = heatmap!(ax3, xc*CharUnits.length.val, zc[1:end-32]*CharUnits.length.val, data[:,1:end-32], colormap = :viridis)
Colorbar(fig[2, 2], hm3, label = "Stress [MPa]")
#sensitivities
ax1 = Axis(fig[1,1], xlabel = "x [km]", ylabel = "z [km]", title = "sensitivites ρ ")
hm1 = heatmap!(ax1, xc*CharUnits.length.val, zc[1:end-18]*CharUnits.length.val, sens_ρ[:,1:end-18], colormap = :ocean, colorrange = (-1.1*maximum(abs.(sens_ρ[:,1:end-32])), 1.5*maximum(abs.(sens_ρ[:,1:end-32]))))
Colorbar(fig[1,2], hm1, label = "sensitivites")

ax4 = Axis(fig[2,1], xlabel = "x [km]", ylabel = "z [km]", title = "sensitivites η")
hm4 = heatmap!(ax4, xc*CharUnits.length.val, zc[1:end-18]*CharUnits.length.val, sens_η[:,1:end-18], colorrange = (-0.5*maximum(abs.(sens_η[:,1:end-32])), 1.2*maximum(abs.(sens_η[:,1:end-32]))), colormap = :ocean)
Colorbar(fig[2, 2], hm4, label = "sensitivites")

ax5 = Axis(fig[3,3], xlabel = "x [km]", ylabel = "z [km]", title = "sensitivities ρ + η")
hm5 = heatmap!(ax5, xc*CharUnits.length.val, zc*CharUnits.length.val[1:end-18], sensitivities[1][:,1:end-18] + average_2D(sensitivities[2])[:,1:end-18], colormap= :ocean, colorrange = (-maximum(abs.(sensitivities[1] + average_2D(sensitivities[2]))), maximum(abs.(sensitivities[1] + average_2D(sensitivities[2])))))
lines!(ax5, [xc[IDX[1][1][1]]*CharUnits.length.val, xc[IDX[end][1][1]]*CharUnits.length.val], [zc[IDX[1][1][2]]*CharUnits.length.val, zc[IDX[end][1][2]]*CharUnits.length.val], color = :red, linewidth=5)
Colorbar(fig[3, 4], hm5, label = "sensitivites")
display(fig)