# example_2.jl
#
# We solve the 2D elliptic PDE:
#   ∂( k * ∂u/∂x )/∂x + ∂( k * ∂u/∂z )/∂z = 0
# with boundary conditions:
#   u(0) = 0, u(L) = 100
#
# This example is a bit more involved and uses the add-on routines, to allow passing 1 or more fields in 2D to the residual routine
# We also use the `nonlinear_solution` routine to solve the system of equations (which converges in 1 step if it is linear)

using SparseDiffTools, LinearAlgebra
using Symbolics
using GLMakie
Makie.inline!(true)

include("addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

"""
    R!(F::AbstractArray, u::AbstractArray, Δ, N, BC)

Compute the residual of the PDE; in this case: 
    R = ∂( k*∂u/∂x )/∂x + ∂( k*∂u/∂z )/∂z 
    with u[1,:] .= u0; u[N,:] .= uEnd
"""
function Res!(Fvec::AbstractVector{T}, UU::Vector{<:AbstractArray{T}}, Δ::NTuple, N::NTuple, BC::NamedTuple, Params::NamedTuple) where T<:Number

    dx, dz = Δ
    
    U = UU[1]       # U is assumed to be the first 
    F = zero(U)

    # ∂( kx*∂u/∂x )/∂x + ∂( kz*∂u/∂z )/∂z = 0
    Nx = N[1][1]
    Nz = N[1][2]
    
    kx = ones(Nx-1,Nz  )
    kz = ones(Nx  ,Nz-1)
    
    # central
    F[2:Nx-1,2:Nz-1] =  diff( kx[:,2:Nz-1].* diff(U[:,2:Nz-1],dims=1)./dx, dims=1 )./dx - 
                        diff( kz[2:Nx-1,:].* diff(U[2:Nx-1,:],dims=2)./dz, dims=2 )./dz 
    # bot & top
    F[:,   1] = U[:,1 ] .- BC.u0
    F[:,  Nz] = U[:,Nz] .- BC.uEnd
    
    # zero flux left and right
    F[1,  2:Nz-1] = (U[2, 2:Nz-1] .- U[1,   2:Nz - 1])/dx .- 0.0
    F[Nx, 2:Nz-1] = (U[Nx,2:Nz-1] .- U[Nx-1,2:Nz - 1])/dx .- 0.0

    Fvec .= F[:]
    return nothing
end

function Res!(Fup::AbstractVector{T}, up::AbstractVector{T}, Δ, N, BC, Params) where {T<:Number}
    return Res!(Fup, vec_2_vecarray(up, N), Δ, N, BC, Params)
end

Res_closed! = (F,U) -> Res!(F, U, Δ, N, BC, Params)        # create a function with only 1 input parameter

# Setup problem
N = ((11,12),)
U = [rand(N[1]...)]        # 
x=range(0,2,N[1][1])
z=range(0,2,N[1][2])

dx = x.step.hi
dz = z.step.hi
Δ  = (dx,dz)
Params      = NamedTuple()                 # additional parameters of the PDE

u0 = 0;
uEnd= 100;
BC = (;u0, uEnd)
U[1][:,1]   .=BC.u0;
U[1][:,end] .=BC.uEnd;
F= vecarray_2_vec(U)

# Sparsity pattern of jacobian
sparsity    =   jacobian_sparsity(Res_closed!,F,U)
J           =   Float64.(sparsity)
colors      =   matrix_colors(J) 

# Solve
Usol = nonlinear_solution(F, U, J, colors);


# plot
heatmap(x,z,Usol[1])