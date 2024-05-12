# example_3.jl
#
# We solve the two elliptic PDEs:
#   ∂( ku * ∂u/∂x )/∂x = 0
#   ∂( kp * ∂p/∂x )/∂x = 0
# with boundary conditions:
#   u(0) = 0, u(L) = 100
#   p(0) = 11, p(L) = 21

using LinearAlgebra, SparseDiffTools
using Symbolics
using GLMakie
Makie.inline!(true)

include("addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

"""
    Res!(Fup::AbstractVector{T}, up::Vector{<:AbstractArray{T}}, Δ, N, BC, Params)

Compute the residual of one or more (coupled) PDEs. `up` is a vector of abstract arrays that contains the values of each of the fields, 
`Δ` is the grid spacing, `N` is the number of grid points, `BC` contain the boundary conditions, and `Params` can contain additional 
parameters of the PDE.
`Fup` is a vector that contains the residual 
"""
function Res!(Fup::AbstractVector{T}, up::Vector{<:AbstractArray{T}}, Δ::NTuple, N::NTuple, BC::NamedTuple, Params::NamedTuple) where T<:Number
    u,p  = up[1], up[2]
    Nu, Np = N[1], N[2] 
    dx = Δ

    ind = cumsum(prod.(N))
    Fu = Fup[       1:ind[1]]
    Fp = Fup[ind[1]+1:ind[2]]

    # residual function of the PDE:
    # ∂( k*∂u/∂x )/∂x = 0
    ku = ones(Nu.-1)
    #k = 1.0 .+ 0.1*u[2:end]./100
    #k = range(1,10,length(u)-1)
    Fu[1]       = u[1] - BC.u0
    Fu[2:Nu-1]  = diff( ku.* diff(u)./dx )./dx; # - 10*sin.(x[2:end-1]) - x[2:end-1].^2.0
    Fu[Nu]      = u[Nu] - BC.uEnd
    
    # residual function of the PDE:
    # ∂( k*∂u/∂x )/∂x = 0
    kp = ones(Np.-1)
    kp = 1.0 .+ 0.1*p[2:end]./100
    #k = range(1,10,length(u)-1)
    Fp[1]       = p[1] - BC.p0
    Fp[2:Np-1]  = diff( kp.* diff(p)./dx )./dx; # - 10*sin.(x[2:end-1]) - x[2:end-1].^2.0
    Fp[Np]      = p[Np] - BC.pEnd
    
    Fup[1:length(Fu[:])]         = Fu[:]
    len                          = length(Fu[:]);
    Fup[len+1:len+length(Fp[:])] = Fp[:]
    
    return nothing
end

function Res!(Fup::AbstractVector{T}, up::AbstractVector{T}, Δ, N, BC, Params) where {T<:Number}
    return Res!(Fup, vec_2_vecarray(up, N), Δ, N, BC, Params)
end
Res_closed! = (F,u) -> Res!(F,u, Δ, N, BC, Params)            # create a function with only 1 input parameter


# Setup
n = 101
u = rand(n)
p = rand(n)
L = 2;
x=range(0,L,n)
dx = x.step.hi

u0, uEnd = 0, 100;
p0, pEnd = 10, 21;

U           = [u,p]                         # tuple with solution fields
N           = (n,n)                         # number of grid points of each field    
Fup         = zeros(length(u) + length(p))  # Residual vector
BC          = (; u0, uEnd, p0, pEnd)        
Δ           = (dx,)                        # grid spacing
Params      = NamedTuple()                 # additional parameters of the PDE

# Sparsity pattern of jacobian
sparsity    =   jacobian_sparsity(Res_closed!,Fup, U)
J           =   Float64.(sparsity)
colors      =   matrix_colors(J) 

# Solve
Usol = nonlinear_solution(Fup, U, J, colors);

fig, ax = plot(x,Usol[2], axis=(xlabel="x",ylabel="u"))
plot!(ax,x,Usol[1],color=:red)
fig
