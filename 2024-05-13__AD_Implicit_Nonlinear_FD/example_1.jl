# example_1.jl
#
# We solve the elliptic PDE:
#   ∂( k * ∂u/∂x )/∂x = 0
# with boundary conditions:
#   u(0) = 0, u(L) = 100
# and assume that k is either constant or k(x)
#
# This example illustrates how we can use automatic differentiaion and sparsity detection to automatically solve the steady state diffusion equations 

using SparseDiffTools, LinearAlgebra
using Symbolics
using GLMakie
Makie.inline!(true)

"""
    R!(F::AbstractVector, u::AbstractVector, dx, N, u0, uEnd)

Compute the residual of the PDE; in this case: 
    R = ∂( k*∂u/∂x )/∂x 
    with u[0] = u0; u[N]=uEnd
"""
function R!(F::AbstractVector, u::AbstractVector, Δ, N, BC)
    # residual function of the PDE:
    # ∂( k*∂u/∂x )/∂x = 0
    dx = Δ[1]
    
    #k = ones(N.-1)
    k = range(1,10,length(u)-1)     # variable k
    F[1] = u[1] - BC.u0
    F[2:N-1] = diff( k.* diff(u)./dx )./dx; 
    F[N] = u[N] - BC.uEnd
    
    return nothing
end
Res_closed! = (F,u) -> R!(F,u, Δ, N, BC)        # create a function with only 1 input parameter

# Setup
N = 11
u = rand(N)
x = range(0,2,N[1])
Δ = (x.step.hi,)

u0, uEnd = 0, 100;
BC       = (; u0, uEnd)
F        = zero(u)

# Sparsity pattern of jacobian
sparsity    =   jacobian_sparsity(Res_closed!,F, u)
J           =   Float64.(sparsity)
colors      =   matrix_colors(J) 

# compute jacobian in an in-place manner
forwarddiff_color_jacobian!(J, Res_closed!, u, colorvec = colors)

# Compute linear solution:
Res_closed!(F,u)    # compute residual
du = J\-F

u = u+du


#=
# So how does sparsity detection work?
# 1) give every variable a unique name
vars = map(Symbolics.variable, eachindex(u))       

# 2) Evaluate residual function
expr = zero(vars)
Res_closed!(expr, vars)

# 3) Compute nonzero coefficients
sparsity = Symbolics.jacobian_sparsity(expr, vars)
=#