# example_4.jl
#
# We solve the two coupled PDE's which simulate magma migration through viscoelastic rocks in the Earth
#    ∂ϕ/∂t  = -ϕ^m *Pe - De*∂Pe/∂t
# De*∂Pe/∂t = ∂( ϕ^n * (∂Pe/∂z + 1) )/∂z + ϕ^m *Pe
#
# where ϕ= melt fraction and Pe = effective pressure. We use a staggered grid formulation 
# Equations are described in https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/98GL52358; 

using LinearAlgebra, SparseDiffTools
using Symbolics
using GLMakie
Makie.inline!(true)

include("addons_multiphysics_AD.jl")    # bunch of helper files to deal with multiple fields

# helper functions
f(z,zk,λ) = exp.( - (z.-zk).^2.0 ./ λ.^2 )
average(x) = (x[2:end] + x[1:end-1])/2

"""
    Res!(F::AbstractVector, U::Vector{<:AbstractArray}, Δ, N, BC, Params)
"""
function Res!(F::AbstractVector{T}, U::Vector{<:AbstractArray{T}}, Δ::NTuple, N::NTuple, BC::NamedTuple, Params::NamedTuple) where T<:Number
    ϕ, Pe   = U[1], U[2]
    Nϕ      = N[1]
    dz      = Δ
    Fϕ      = zeros(eltype(ϕ),  Nϕ  )
    Fp      = zeros(eltype(Pe), Nϕ+1)
    n       = Params.n
    m       = Params.m   
    ϕold    = Params.ϕold
    Pe_old  = Params.Pe_old
    
    # residual function 1
    # (ϕ-ϕold)/Δt + ϕ^m *Pe + De*(Pe-Pe_old)/Δt 
    Pe_c        = average(Pe)   
    Pe_c_old    = average(Pe_old)   
    Fϕ          = (ϕ - ϕold)./Params.Δt + Params.De*(Pe_c-Pe_c_old)/Δt  + ϕ.^m .*Pe_c
    
    # residual function 2:
    # De*∂Pe/∂t - ∂( ϕ^n * (∂Pe/∂z + 1) )/∂z + ϕ^m *Pe
    Fp[1]       = Pe[1] - BC.p0
    ϕ_c         = average(ϕ)  
    Fp[2:Nϕ]    = Params.De*(Pe[2:end-1] - Pe_old[2:end-1])/Δt -   diff( ϕ.^n .* (diff(Pe)./dz .+ 1))./dz + ϕ_c.^m .* Pe[2:end-1];
    Fp[end]     = Pe[end] - BC.pEnd
    
    F           .= [Fϕ[:]; Fp[:]]
    return nothing
end

function Res!(F::AbstractVector{T}, U::AbstractVector{T}, Δ, N, BC, Params) where {T<:Number}
    return Res!(F, vec_2_vecarray(U, N), Δ, N, BC, Params)
end
Res_closed! = (F,U) -> Res!(F, U, Δ, N, BC, Params)            # create a function with only 1 input parameter

# Setup
Nz =    1001
ϕ  =    rand(Nz)
Pe =    rand(Nz+1)
L  =    2;
z  =    range(-20,150,Nz+1)
zc =    (z[2:end] + z[1:end-1]) ./ 2
dz =    z.step.hi

p0, pEnd = 0, 0;

U           = [ϕ,Pe]                         # tuple with solution fields
N           = (Nz,Nz+1)                        # number of grid points of each field    
F           = zeros(length(ϕ) + length(Pe))  # Residual vector
BC          = (; p0, pEnd)        
Δ           = (dz,)                          # grid spacing

# initial conditions
De          =   1e2;
ϕ0          =   1.0
Δϕ1         =   8.0
Δϕ2         =   1.0
ΔPe1        =   Δϕ1/De
ΔPe2        =   Δϕ2/De
z1          =   0.0
z2          =   40
λ           =   1.0

ϕ           =   ϕ0 .+ Δϕ1*f(zc,z1,λ) .+ Δϕ2*f(zc,z2,λ)
ϕold        =   copy(ϕ)                                  # old melt fraction
Pe          =   -ΔPe1*f(z,z1,λ) .- ΔPe2*f(z,z2,λ)
Pe_old      =   -ΔPe1*f(z,z1,λ) .- ΔPe2*f(z,z2,λ)

m           =   2;
n           =   3; 
Δt          =   1e-2    # good value for De=1e-2
Params      =   (;ϕold,Pe_old, m,n,Δt, De)                # additional parameters of the PDE

# Sparsity pattern of jacobian
sparsity    =   jacobian_sparsity(Res_closed!,F, U)
J           =   Float64.(sparsity)
colors      =   matrix_colors(J) 

function perform_timestepping(Params, U, F, J, colors; tmax=1.0, step_plot=10, Prange=0.03, ϕmax=4, filename= "animation.mp4")
    time = 0.0
    it = 0
    
    fig = Figure(size = (600, 800))
    ax1 = Axis(fig[1, 1], xlabel="melt fraction ϕ", ylabel="Depth Z ", title="time=$time", limits=(0, ϕmax, -10,150))
    ax2 = Axis(fig[1, 2], xlabel="effective pressure Pe", limits=(-Prange,Prange, -10,150))  # create a plot
    
    record(fig, filename) do io
        while time<tmax

            # Solve nonlinear system of equations
            Usol = nonlinear_solution(F, U, J, colors);

            # Update old parameters
            Params.ϕold     .= copy(Usol[1])
            Params.Pe_old   .= copy(Usol[2])
            U               .= Usol;

            time += Δt
            it += 1
            println("it $it, time=$time")

            # Create a plot
            if mod(it,step_plot) == 0
                empty!(ax1)
                empty!(ax2)
                lines!(ax2, Usol[2], z,  color=:blue)
                lines!(ax1, Usol[1], zc, color=:black)

                #
                ax1.title=rpad("De=$De, time=$(round(time, digits=3))",20)
                recordframe!(io)  # record a new frame
            end

        end
    end

    return U

end

# De=1e-2 (viscous) (note: plotting requires different ranges)
#Usol = perform_timestepping(Params, U, F, J, colors, tmax=25, step_plot=10, Prange=1.5, ϕmax=8, filename="De1e_2_1001points.mp4")

#De=1e2 (elastic)
Usol = perform_timestepping(Params, U, F, J, colors, tmax=25, step_plot=10, Prange=0.03, ϕmax=2.5, filename="De1e2_1001points.mp4")