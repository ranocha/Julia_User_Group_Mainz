# example_2.jl
#
# We solve the 2D stokes equations using a staggered grid formulation with Vx,Vz,P as unknowns
# The equations are:
#        ∂Vx/∂x + ∂Vz/∂z = -P/K
#   -∂P/∂x + ∂τxx/∂x + ∂τxz/∂x = 0
#   -∂P/∂z + ∂τxz/∂x + ∂τzz/∂x = ρg
#
#   τxx = 2 η εxx 
#   τzz = 2 η εzz 
#   τxz =   η εxz 
#
#   εxx = ∂Vx/∂x
#   εzz = ∂Vz/∂z
#   εxz = (∂Vx/∂z + ∂Vz/∂x)/2

using SparseDiffTools, LinearAlgebra
using Symbolics
using GLMakie
Makie.inline!(true)

include("addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

"""
    Res!(F::AbstractArray, u::AbstractArray, Δ, N, BC)
"""
function Res!(Fvec::AbstractVector{T}, U::Vector{<:AbstractArray{T}}, Δ::NTuple, N::NTuple, BC::NamedTuple, Params::NamedTuple) where T<:Number

    dx, dz      = Δ         # grid spacing
    Vx, Vz, P   = U         # extract 2D fields
    Nx, Nz      = N[3]      # grid size

    # Strainrates
    εxx = diff(Vx,dims=1)./dx
    εzz = diff(Vz,dims=2)./dz
    εxz = zeros(eltype(Vx), Nx+1, Nz+1)
    εxz[2:end-1,2:end-1] = (diff(Vx[2:end-1,:],dims=2)./dz + diff(Vz[:,2:end-1],dims=1)./dx)

    εxz2_c = (εxz[2:end,2:end].^2 + εxz[1:end-1,2:end].^2 + εxz[1:end-1,1:end-1].^2 + εxz[2:end,1:end-1].^2)/4
    εII_c = sqrt.(0.5*εxx.^2 + 0.5.*εzz.^2 + εxz2_c.^2) # 2nd invariant @ center

    # Update deviatoric stress
    τxx = 2*Params.ηc.*εxx
    τzz = 2*Params.ηc.*εzz
    τxz = 2*Params.ηv.*εxz

    # Free slip BC's
    τxz[1,  :] .= 0
    τxz[:,end] .= 0
    τxz[:,  1] .= 0
    τxz[end,:] .= 0
    
    # Mass balance
    Fmass = zero(P)
    Fmass = εxx +  εzz + K.*P 

    # x-momentum balance
    Fforce1             = zero(Vx)
    Fforce1[2:end-1,:]  = diff(-P + τxx, dims=1)./dx + diff(τxz[2:end-1,:], dims=2)./dz
    Fforce1[1,:]        = Vx[1,:] .- BC.εbg*Params.x[1];
    Fforce1[end,:]      = Vx[end,:] .- BC.εbg*Params.x[end];

    # z-momentum balance
    Fforce2             =   zero(Vz)
    Fforce2[:,2:end-1]  =   diff(τxz[:,2:end-1], dims=1)./dx + diff(-P + τzz, dims=2)./dz - Params.g*(Params.ρ[:,2:end] + Params.ρ[:,1:end-1])/2 
    Fforce2[:,      1]  =   Vz[:,1] .+ BC.εbg*Params.z[1];
    Fforce2[:,    end]  =   Vz[:,end] .+ BC.εbg*Params.z[end];

    Fvec .= [Fforce1[:]; Fforce2[:]; Fmass[:]]
    return nothing
end

function Res!(Fup::AbstractVector{T}, up::AbstractVector{T}, Δ, N, BC, Params) where {T<:Number}
    return Res!(Fup, vec_2_vecarray(up, N), Δ, N, BC, Params)
end
Res_closed! = (F,U) -> Res!(F, U, Δ, N, BC, Params)        # create a function with only 1 input parameter

# Setup problem
Nx, Nz = 201,201

N   =   ((Nx+1,Nz), (Nx,Nz+1), (Nx,Nz))
Vx  =   rand(N[1]...)
Vz  =   rand(N[2]...)
P   =   rand(N[3]...)
x   =   range(-1,1,Nx+1)
z   =   range(-1,1,Nz+1)
xc  =   (x[2:end] .+ x[1:end-1])/2
zc  =   (z[2:end] .+ z[1:end-1])/2

# Material fields
ρ   =    ones(Nx,Nz)             # density
ηc  =    ones(Nx,Nz)             # viscosity @ center
ηv  =    ones(Nx+1,Nz+1)    # viscosity @ center
K   =    1e-10
g   =    1;
εbg =    0;

# Helper fields
εxx =  zeros(Nx,Nz)
εzz =  zeros(Nx,Nz)
εxz =  zeros(Nx+1,Nz+1)
τxx =  zeros(Nx,Nz)
τzz =  zeros(Nx,Nz)
τxz =  zeros(Nx+1,Nz+1)

dx  = x.step.hi
dz  = z.step.hi
Δ   = (dx,dz)

# Set density/viscosity anomaly
for (i,xv) in enumerate(xc), (j,zv) in enumerate(zc)
    R = 0.3
    if (xv^2 + zv^2) < R^2
        ρ[i,j]  = 2.0
        ηc[i,j] = 2.0
    end
end


Params  = (; K, εxx,εzz,εxz,τxx,τzz,τxz, ρ, ηc, ηv, g,xc,zc,x,z)                 # additional parameters of the PDE
BC      = (; εbg)
U       = [Vx,Vz,P]
F       = vecarray_2_vec(U)


# Sparsity pattern of jacobian
@time sparsity    =   jacobian_sparsity(Res_closed!,F,U)
J           =   Float64.(sparsity)
colors      =   matrix_colors(J) 


# Solve system of equations:
@time Usol = nonlinear_solution(F, U, J, colors, maxit=1);

# plot
fig, ax, hm = heatmap(x,z,Usol[2])
Colorbar(fig[1,2], limits = extrema(Usol[2]), label="Vz")
fig
